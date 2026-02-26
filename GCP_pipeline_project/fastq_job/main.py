import os
import sys
import json
import uuid
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

from google.cloud import storage
from google.cloud import bigquery

# -----------------------------
# Helpers: GCS
# -----------------------------
def parse_gs_uri(uri: str) -> Tuple[str, str]:
    if not uri.startswith("gs://"):
        raise ValueError(f"Not a gs:// URI: {uri}")
    parts = uri[5:].split("/", 1)
    bucket = parts[0]
    blob = parts[1] if len(parts) > 1 else ""
    return bucket, blob

def download_blob(gs_uri: str, local_path: str) -> None:
    bucket_name, blob_name = parse_gs_uri(gs_uri)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    Path(local_path).parent.mkdir(parents=True, exist_ok=True)
    blob.download_to_filename(local_path)

def upload_text(bucket_name: str, blob_name: str, text: str, content_type: str) -> None:
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.upload_from_string(text, content_type=content_type)

# -----------------------------
# Helpers: BigQuery
# -----------------------------
def bq_insert_rows(project: str, dataset: str, table: str, rows: List[dict]) -> None:
    if not rows:
        return
    client = bigquery.Client(project=project)
    table_id = f"{project}.{dataset}.{table}"
    errors = client.insert_rows_json(table_id, rows)
    if errors:
        raise RuntimeError(f"BigQuery insert errors into {table_id}: {errors}")

# -----------------------------
# Bio/QC (minimal for now)
# -----------------------------
def fastq_basic_qc(fastq_path: str) -> Dict:
    # very small QC: counts reads and total bases assuming 4-line FASTQ records
    reads = 0
    total_bases = 0
    min_len = None
    max_len = 0

    with open(fastq_path, "r", encoding="utf-8", errors="ignore") as f:
        i = 0
        for line in f:
            i += 1
            if i % 4 == 2:  # sequence line
                seq = line.strip()
                L = len(seq)
                reads += 1
                total_bases += L
                min_len = L if min_len is None else min(min_len, L)
                max_len = max(max_len, L)

    mean_len = (total_bases / reads) if reads else 0.0
    return {
        "reads": reads,
        "total_bases": total_bases,
        "read_length": {"min": min_len or 0, "mean": mean_len, "max": max_len},
    }

# -----------------------------
# BWA + samtools + featureCounts
# -----------------------------
def run(cmd: List[str], stdout_path: str | None = None) -> str:
    import subprocess

    try:
        if stdout_path:
            with open(stdout_path, "w", encoding="utf-8") as out:
                p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=False)
        else:
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)

        if p.returncode != 0:
            raise RuntimeError(
                "Command failed:\n"
                f"  {' '.join(cmd)}\n"
                f"STDERR:\n{p.stderr}\n"
                f"STDOUT:\n{p.stdout if p.stdout else ''}"
            )

        return p.stdout if p.stdout else ""
    except Exception as e:
        raise

def download_reference_and_index(ref_uri: str, workdir: str) -> str:
    """
    Downloads reference.fasta and the BWA index sidecars if present.
    Returns local path to reference fasta.
    """
    Path(workdir).mkdir(parents=True, exist_ok=True)
    local_ref = str(Path(workdir) / "reference.fasta")
    download_blob(ref_uri, local_ref)

    # try to download BWA index sidecars if they exist next to the fasta
    # (your bucket already has .amb/.ann/.bwt/.pac/.sa)
    suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    for suf in suffixes:
        try:
            download_blob(ref_uri + suf, local_ref + suf)
        except Exception:
            # If missing, BWA can still build indexes (slow), but for toy runs it's fine.
            pass

    return local_ref

def bwa_align_sorted_bam(fastq_path: str, ref_fasta: str, out_prefix: str) -> Tuple[str, str]:
    sam = out_prefix + ".sam"
    bam = out_prefix + ".bam"
    sorted_bam = out_prefix + ".sorted.bam"
    idxstats = out_prefix + ".idxstats.tsv"

    # bwa mem -> stdout, so write stdout to sam
    run(["bwa", "mem", "-t", "2", ref_fasta, fastq_path], stdout_path=sam)

    run(["samtools", "view", "-bS", sam, "-o", bam])
    run(["samtools", "sort", "-o", sorted_bam, bam])
    run(["samtools", "index", sorted_bam])

    # idxstats -> stdout, so write stdout to idxstats file
    run(["samtools", "idxstats", sorted_bam], stdout_path=idxstats)

    return sorted_bam, idxstats

def parse_idxstats(idxstats_path: str) -> List[dict]:
    rows = []
    with open(idxstats_path, "r", encoding="utf-8") as f:
        for line in f:
            contig, length, mapped, unmapped = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "contig": contig,
                    "mapped": int(mapped),
                    "unmapped": int(unmapped),
                }
            )
    return rows

def run_featurecounts(gtf_path: str, bam_path: str, out_tsv: str) -> None:
    # featureCounts writes a header block; still fine, we parse later
    run([
        "featureCounts",
        "-a", gtf_path,
        "-o", out_tsv,
        "-t", "exon",
        "-g", "gene_id",
        bam_path
    ])

def parse_featurecounts_tsv(out_tsv: str) -> List[Tuple[str, int]]:
    """
    Returns list of (gene_id, count). Handles typical featureCounts output.
    The last column is the BAM count.
    """
    results: List[Tuple[str, int]] = []
    with open(out_tsv, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            if line.startswith("Geneid"):
                continue
            parts = line.rstrip("\n").split("\t")
            # Standard columns: Geneid Chr Start End Strand Length <bam>
            gene_id = parts[0]
            count = int(parts[-1])
            results.append((gene_id, count))
    return results

# -----------------------------
# Main
# -----------------------------
def main() -> None:
    input_uri = os.environ.get("INPUT_URI") or (sys.argv[1] if len(sys.argv) > 1 else None)
    results_bucket = os.environ.get("RESULTS_BUCKET")

    ref_uri = os.environ.get("REF_FASTA_URI")
    gtf_uri = os.environ.get("GTF_URI")

    if not input_uri or not results_bucket:
        raise SystemExit("Missing INPUT_URI and/or RESULTS_BUCKET")

    if not ref_uri:
        raise SystemExit("Missing REF_FASTA_URI (gs://... reference fasta)")
    if not gtf_uri:
        raise SystemExit("Missing GTF_URI (gs://... annotation.gtf)")

    # BigQuery config
    bq_project = os.environ.get("BQ_PROJECT", "fastq-data-pipeline")
    bq_dataset = os.environ.get("BQ_DATASET", "fastq_qc")
    qc_table = os.environ.get("BQ_TABLE", "qc_runs")  # keep your existing qc table name

    # derive sample name from object name
    _, object_name = parse_gs_uri(input_uri)
    sample_name = object_name.rsplit("/", 1)[-1]

    run_id = str(uuid.uuid4())
    ran_at = datetime.now(timezone.utc)

    # workspace
    workdir = "/tmp/work"
    Path(workdir).mkdir(parents=True, exist_ok=True)

    local_fastq = str(Path(workdir) / "input.fastq")
    local_gtf = str(Path(workdir) / "annotation.gtf")

    download_blob(input_uri, local_fastq)
    download_blob(gtf_uri, local_gtf)

    local_ref = download_reference_and_index(ref_uri, workdir)
    sorted_bam, idxstats_path = bwa_align_sorted_bam(local_fastq, local_ref, str(Path(workdir) / "sample"))

    # Upload idxstats + contig counts json
    upload_text(
        results_bucket,
        f"{sample_name}.idxstats.tsv",
        Path(idxstats_path).read_text(),
        "text/tab-separated-values",
    )

    contig_counts = parse_idxstats(idxstats_path)
    upload_text(
        results_bucket,
        f"{sample_name}.counts.json",
        json.dumps(contig_counts, indent=2),
        "application/json",
    )

    # featureCounts
    gene_tsv = str(Path(workdir) / f"{sample_name}.gene_counts.tsv")
    run_featurecounts(local_gtf, sorted_bam, gene_tsv)

    upload_text(
        results_bucket,
        f"{sample_name}.gene_counts.tsv",
        Path(gene_tsv).read_text(),
        "text/tab-separated-values",
    )

    gene_counts = parse_featurecounts_tsv(gene_tsv)
    upload_text(
        results_bucket,
        f"{sample_name}.gene_counts.json",
        json.dumps([{ "gene_id": g, "count": c } for g, c in gene_counts], indent=2),
        "application/json",
    )

    # QC summary (toy/basic)
    qc = fastq_basic_qc(local_fastq)

    # Insert qc row
    qc_row = {
        "run_id": run_id,
        "ran_at_utc": ran_at.isoformat(),
        "input_uri": input_uri,
        "sample_name": sample_name,
        "reads": int(qc["reads"]),
        "total_bases": int(qc["total_bases"]),
        "read_len_min": int(qc["read_length"]["min"]),
        "read_len_mean": float(qc["read_length"]["mean"]),
        "read_len_max": int(qc["read_length"]["max"]),
        # keep any extra columns in your qc_runs schema aligned as needed
    }
    bq_insert_rows(bq_project, bq_dataset, qc_table, [qc_row])

    # Insert contig counts (fastq_qc.contig_counts)
    contig_rows = []
    for r in contig_counts:
        contig_rows.append({
            "run_id": run_id,
            "ran_at_utc": ran_at.isoformat(),
            "input_uri": input_uri,
            "sample_name": sample_name,
            "contig": r["contig"],
            "mapped": r["mapped"],
            "unmapped": r["unmapped"],
        })
    bq_insert_rows(bq_project, bq_dataset, "contig_counts", contig_rows)

    # Insert gene counts (fastq_qc.gene_counts)
    gene_rows = []
    for gene_id, count in gene_counts:
        gene_rows.append({
            "run_id": run_id,
            "ran_at_utc": ran_at.isoformat(),
            "input_uri": input_uri,
            "sample_name": sample_name,
            "gene_id": gene_id,
            "count": int(count),
        })
    bq_insert_rows(bq_project, bq_dataset, "gene_counts", gene_rows)

    print(f"OK run_id={run_id} sample_name={sample_name}")
    print(f"Wrote gs://{results_bucket}/{sample_name}.gene_counts.tsv")
    print(f"Inserted BigQuery rows into {bq_project}.{bq_dataset}.gene_counts")

if __name__ == "__main__":
    main()
