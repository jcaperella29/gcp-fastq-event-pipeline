import csv
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

from google.cloud import storage

from db import get_db_connection, mark_run_complete, mark_run_failed, mark_run_running

BARCODE_TO_SAMPLE = {
    "ACGTACGT": "sample1",
    "TGCATGCA": "sample2",
    "GATTACAG": "sample3",
}
BARCODE_LENGTH = 8


def parse_gcs_uri(uri: str):
    parsed = urlparse(uri)
    if parsed.scheme != "gs":
        raise ValueError(f"Invalid GCS URI: {uri}")
    return parsed.netloc, parsed.path.lstrip("/")


def run_cmd(cmd: list[str], cwd: str | None = None):
    print("RUN:", " ".join(cmd))
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr)
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(cmd)}")
    return result


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gtf_gene_intervals(gtf_path: str):
    gene_intervals = defaultdict(list)

    with open_text(gtf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue

            contig, _, feature, start, end, _, strand, _, attrs = fields
            if feature != "exon":
                continue

            gene_id = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id "):
                    gene_id = attr.split('"')[1]
                    break

            if gene_id is None:
                continue

            gene_intervals[gene_id].append(
                {
                    "chrom": contig,
                    "start": int(start) - 1,
                    "end": int(end),
                    "strand": strand,
                }
            )

    return dict(gene_intervals)


def compute_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def list_fastq_uris(data_dir: str):
    if data_dir.startswith("gs://"):
        storage_client = storage.Client()
        bucket_name, prefix = parse_gcs_uri(data_dir)
        bucket = storage_client.bucket(bucket_name)
        blobs = list(bucket.list_blobs(prefix=prefix))

        return sorted(
            [
                f"gs://{bucket_name}/{b.name}"
                for b in blobs
                if b.name.endswith(".fastq")
                or b.name.endswith(".fq")
                or b.name.endswith(".fastq.gz")
                or b.name.endswith(".fq.gz")
            ]
        )

    return sorted(
        [
            str(p)
            for p in Path(data_dir).iterdir()
            if p.name.endswith(".fastq")
            or p.name.endswith(".fq")
            or p.name.endswith(".fastq.gz")
            or p.name.endswith(".fq.gz")
        ]
    )


def download_gcs_file(uri: str, local_path: str):
    storage_client = storage.Client()
    bucket_name, blob_name = parse_gcs_uri(uri)
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.download_to_filename(local_path)


def stage_inputs_locally(
    reference_fasta: str,
    gtf_file: str,
    fastq_uris: list[str],
    work_dir: str,
):
    ref_local = os.path.join(work_dir, "reference.fasta")
    gtf_local = os.path.join(work_dir, "annotation.gtf")
    fastq_dir = os.path.join(work_dir, "raw_fastqs")
    os.makedirs(fastq_dir, exist_ok=True)

    if reference_fasta.startswith("gs://"):
        download_gcs_file(reference_fasta, ref_local)
    else:
        shutil.copy(reference_fasta, ref_local)

    if gtf_file.startswith("gs://"):
        download_gcs_file(gtf_file, gtf_local)
    else:
        shutil.copy(gtf_file, gtf_local)

    fastq_locals = []
    for i, fq_uri in enumerate(fastq_uris, start=1):
        suffix = ".fastq.gz" if fq_uri.endswith(".gz") else ".fastq"
        local_path = os.path.join(fastq_dir, f"input_{i}{suffix}")
        if fq_uri.startswith("gs://"):
            download_gcs_file(fq_uri, local_path)
        else:
            shutil.copy(fq_uri, local_path)
        fastq_locals.append(local_path)

    return ref_local, gtf_local, fastq_locals


def open_fastq(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def demultiplex_fastqs(local_fastqs: list[str], out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    sample_fastq_paths = {
        sample: os.path.join(out_dir, f"{sample}.fastq")
        for sample in sorted(set(BARCODE_TO_SAMPLE.values()))
    }

    handles = {sample: open(path, "w") for sample, path in sample_fastq_paths.items()}

    metrics = {
        "total_reads": 0,
        "total_bases": 0,
        "gc_bases": 0,
        "barcode_unassigned_reads": 0,
        "reads_per_sample": {sample: 0 for sample in sample_fastq_paths},
    }

    try:
        for fastq_path in local_fastqs:
            with open_fastq(fastq_path) as handle:
                while True:
                    header = handle.readline()
                    if not header:
                        break

                    seq = handle.readline().strip()
                    plus = handle.readline()
                    qual = handle.readline().strip()

                    if not seq or not plus or not qual:
                        raise ValueError(f"Malformed FASTQ record encountered in {fastq_path}")

                    metrics["total_reads"] += 1
                    metrics["total_bases"] += len(seq)
                    metrics["gc_bases"] += sum(1 for b in seq.upper() if b in {"G", "C"})

                    if len(seq) <= BARCODE_LENGTH:
                        metrics["barcode_unassigned_reads"] += 1
                        continue

                    barcode = seq[:BARCODE_LENGTH]
                    trimmed_seq = seq[BARCODE_LENGTH:]
                    trimmed_qual = qual[BARCODE_LENGTH:]

                    sample_name = BARCODE_TO_SAMPLE.get(barcode)
                    if sample_name is None:
                        metrics["barcode_unassigned_reads"] += 1
                        continue

                    out = handles[sample_name]
                    out.write(header)
                    out.write(trimmed_seq + "\n")
                    out.write(plus)
                    out.write(trimmed_qual + "\n")
                    metrics["reads_per_sample"][sample_name] += 1
    finally:
        for h in handles.values():
            h.close()

    avg_read_length = (
        metrics["total_bases"] / metrics["total_reads"] if metrics["total_reads"] else 0.0
    )
    gc_percent = (
        100.0 * metrics["gc_bases"] / metrics["total_bases"] if metrics["total_bases"] else 0.0
    )

    metrics["avg_read_length"] = round(avg_read_length, 2)
    metrics["gc_percent"] = round(gc_percent, 2)

    return sample_fastq_paths, metrics


def align_and_sort_bam(sample_fastq_paths: dict[str, str], reference_fasta: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    run_cmd(["bwa", "index", reference_fasta])

    bam_paths = {}
    for sample, fastq_path in sample_fastq_paths.items():
        if os.path.getsize(fastq_path) == 0:
            continue

        sam_path = os.path.join(out_dir, f"{sample}.sam")
        bam_path = os.path.join(out_dir, f"{sample}.sorted.bam")

        with open(sam_path, "w") as sam_fh:
            result = subprocess.run(
                ["bwa", "mem", reference_fasta, fastq_path],
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(f"bwa mem failed for {sample}: {result.stderr}")
            sam_fh.write(result.stdout)

        run_cmd(["samtools", "sort", "-o", bam_path, sam_path])
        run_cmd(["samtools", "index", bam_path])

        bam_paths[sample] = bam_path

    return bam_paths


def run_featurecounts(gtf_file: str, bam_paths: dict[str, str], out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    counts_txt = os.path.join(out_dir, "featurecounts.txt")

    bam_list = [bam_paths[s] for s in sorted(bam_paths)]
    if not bam_list:
        raise ValueError("No BAM files were produced for featureCounts.")

    cmd = [
        "featureCounts",
        "-a",
        gtf_file,
        "-o",
        counts_txt,
        "-g",
        "gene_id",
        "-t",
        "exon",
    ] + bam_list

    run_cmd(cmd)
    return counts_txt


def parse_featurecounts_matrix(counts_txt: str):
    with open(counts_txt, "r") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]

    data_lines = [line for line in lines if not line.startswith("#")]
    if len(data_lines) < 2:
        raise ValueError("featureCounts output is missing expected data.")

    header = data_lines[0].split("\t")
    rows = [line.split("\t") for line in data_lines[1:]]

    bam_columns = header[6:]
    sample_order = []
    bam_to_sample = {}

    for bam_path in bam_columns:
        bam_name = Path(bam_path).name
        sample = bam_name.replace(".sorted.bam", "").replace(".bam", "")
        sample_order.append(sample)
        bam_to_sample[bam_path] = sample

    gene_names = []
    count_matrix = {sample: {} for sample in sample_order}

    for row in rows:
        gene_id = row[0]
        gene_names.append(gene_id)

        for bam_path, count_str in zip(bam_columns, row[6:]):
            sample = bam_to_sample[bam_path]
            count_matrix[sample][gene_id] = int(float(count_str))

    return count_matrix, sample_order, gene_names


def build_bed_from_bams(bam_paths: dict[str, str], gene_intervals: dict, bed_path: str):
    with open(bed_path, "w") as out:
        for sample, bam_path in bam_paths.items():
            result = subprocess.run(
                ["samtools", "view", bam_path],
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(f"samtools view failed for {bam_path}: {result.stderr}")

            for line in result.stdout.splitlines():
                fields = line.split("\t")
                if len(fields) < 11:
                    continue

                read_name = fields[0]
                flag = int(fields[1])
                chrom = fields[2]
                pos_1based = int(fields[3])
                seq = fields[9]

                if chrom == "*" or seq == "*":
                    continue

                start = pos_1based - 1
                end = start + len(seq)
                strand = "-" if (flag & 16) else "+"

                best_gene = None
                best_overlap = 0
                for gene_id, intervals in gene_intervals.items():
                    for interval in intervals:
                        if interval["chrom"] != chrom:
                            continue
                        ov = compute_overlap(start, end, interval["start"], interval["end"])
                        if ov > best_overlap:
                            best_overlap = ov
                            best_gene = gene_id

                name = f"{sample}|{best_gene or 'unassigned'}|{read_name}"
                out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")


def write_count_matrix_csv(count_matrix: dict, sample_names: list[str], gene_names: list[str], output_dir: str, run_id: int):
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"count_matrix_run_{run_id}.csv")

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample"] + gene_names)
        for sample in sample_names:
            writer.writerow([sample] + [count_matrix.get(sample, {}).get(gene, 0) for gene in gene_names])

    return out_path


def upload_file_to_gcs(local_path: str, bucket_name: str, object_name: str, content_type: str | None = None):
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(object_name)
    if content_type:
        blob.upload_from_filename(local_path, content_type=content_type)
    else:
        blob.upload_from_filename(local_path)


def main():
    run_id_raw = os.getenv("RUN_ID")
    if not run_id_raw:
        raise RuntimeError("RUN_ID environment variable is missing")
    run_id = int(run_id_raw)

    results_bucket = os.environ["RESULTS_BUCKET"]
    data_dir = os.environ.get("DATA_DIR", "../data")
    reference_fasta = os.environ.get("REFERENCE_FASTA", "../reference.fasta")
    gtf_file = os.environ.get("ANNOTATION_GTF", "../annotation.gtf")

    conn = get_db_connection()

    try:
        mark_run_running(conn, run_id)

        fastq_uris = list_fastq_uris(data_dir)
        if not fastq_uris:
            raise FileNotFoundError(f"No FASTQ files found in {data_dir}")

        with tempfile.TemporaryDirectory() as work_dir:
            ref_local, gtf_local, fastq_locals = stage_inputs_locally(
                reference_fasta=reference_fasta,
                gtf_file=gtf_file,
                fastq_uris=fastq_uris,
                work_dir=work_dir,
            )

            demux_dir = os.path.join(work_dir, "demux")
            align_dir = os.path.join(work_dir, "alignments")
            counts_dir = os.path.join(work_dir, "counts")
            os.makedirs(counts_dir, exist_ok=True)

            sample_fastq_paths, demux_metrics = demultiplex_fastqs(fastq_locals, demux_dir)
            bam_paths = align_and_sort_bam(sample_fastq_paths, ref_local, align_dir)
            counts_txt = run_featurecounts(gtf_local, bam_paths, counts_dir)

            count_matrix, sample_names, gene_names = parse_featurecounts_matrix(counts_txt)
            gene_intervals = parse_gtf_gene_intervals(gtf_local)

            bed_local = os.path.join(counts_dir, f"reads_run_{run_id}.bed")
            build_bed_from_bams(bam_paths, gene_intervals, bed_local)

            count_csv_local = write_count_matrix_csv(
                count_matrix=count_matrix,
                sample_names=sample_names,
                gene_names=gene_names,
                output_dir=counts_dir,
                run_id=run_id,
            )

            assigned_reads = sum(
                count_matrix[sample].get(gene, 0)
                for sample in sample_names
                for gene in gene_names
            )

            result = {
                "reads": demux_metrics["total_reads"],
                "assigned_reads": assigned_reads,
                "unassigned_reads": demux_metrics["total_reads"] - assigned_reads,
                "barcode_unassigned_reads": demux_metrics["barcode_unassigned_reads"],
                "total_bases": demux_metrics["total_bases"],
                "avg_read_length": demux_metrics["avg_read_length"],
                "gc_percent": demux_metrics["gc_percent"],
                "data_dir": data_dir,
                "ran_at_utc": datetime.now(timezone.utc).isoformat(),
                "run_id": run_id,
                "fastq_files_processed": fastq_uris,
                "reads_per_sample_after_barcode": demux_metrics["reads_per_sample"],
                "bam_samples": sorted(bam_paths.keys()),
            }

            output_name = f"run_{run_id}.json"

            if results_bucket.startswith("gs://"):
                results_bucket_name = results_bucket.replace("gs://", "", 1)

                qc_local = os.path.join(counts_dir, output_name)
                with open(qc_local, "w") as fh:
                    json.dump(result, fh, indent=2)

                upload_file_to_gcs(
                    qc_local,
                    results_bucket_name,
                    f"qc-results/{output_name}",
                    content_type="application/json",
                )
                upload_file_to_gcs(
                    count_csv_local,
                    results_bucket_name,
                    f"count-matrices/count_matrix_run_{run_id}.csv",
                )
                upload_file_to_gcs(
                    bed_local,
                    results_bucket_name,
                    f"beds/reads_run_{run_id}.bed",
                )

                output_uri = f"gs://{results_bucket_name}/qc-results/{output_name}"
            else:
                os.makedirs(results_bucket, exist_ok=True)

                qc_local = os.path.join(results_bucket, output_name)
                with open(qc_local, "w") as fh:
                    json.dump(result, fh, indent=2)

                shutil.copy(count_csv_local, os.path.join(results_bucket, f"count_matrix_run_{run_id}.csv"))
                shutil.copy(bed_local, os.path.join(results_bucket, f"reads_run_{run_id}.bed"))
                output_uri = qc_local

        print("Processed FASTQs:")
        for fq in fastq_uris:
            print(f"  - {fq}")

        print("Merged count matrix:")
        print(count_matrix)
        print(json.dumps(result, indent=2))

        mark_run_complete(
            conn=conn,
            run_id=run_id,
            read_count=result["reads"],
            total_bases=result["total_bases"],
            avg_read_length=result["avg_read_length"],
            gc_percent=result["gc_percent"],
            output_uri=output_uri,
        )

    except Exception as exc:
        try:
            mark_run_failed(conn, run_id, str(exc))
        except Exception as db_exc:
            print(f"Failed to update DB failure status: {db_exc}")
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
