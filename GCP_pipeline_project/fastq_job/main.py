import gzip
import json
import os
import tempfile
from datetime import datetime, timezone
from urllib.parse import urlparse

from google.cloud import storage

from db import get_db_connection, mark_run_complete, mark_run_failed, mark_run_running


def parse_gcs_uri(uri: str):
    parsed = urlparse(uri)
    if parsed.scheme != "gs":
        raise ValueError(f"Invalid GCS URI: {uri}")
    bucket = parsed.netloc
    blob_name = parsed.path.lstrip("/")
    return bucket, blob_name


def open_fastq(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def compute_fastq_metrics(local_path: str):
    reads = 0
    total_bases = 0
    gc_bases = 0

    with open_fastq(local_path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break

            seq = handle.readline().strip()
            plus = handle.readline()
            qual = handle.readline().strip()

            if not seq or not plus or not qual:
                raise ValueError("Malformed FASTQ record encountered.")

            reads += 1
            total_bases += len(seq)
            gc_bases += sum(1 for base in seq.upper() if base in {"G", "C"})

    avg_read_length = (total_bases / reads) if reads else 0.0
    gc_percent = ((gc_bases / total_bases) * 100.0) if total_bases else 0.0

    return {
        "reads": reads,
        "total_bases": total_bases,
        "avg_read_length": round(avg_read_length, 2),
        "gc_percent": round(gc_percent, 2),
    }


def main():
    input_uri = os.environ["INPUT_URI"]
    run_id = int(os.environ["RUN_ID"])
    results_bucket = os.environ["RESULTS_BUCKET"]

    conn = get_db_connection()
    try:
        mark_run_running(conn, run_id)

        storage_client = storage.Client()
        in_bucket_name, in_blob_name = parse_gcs_uri(input_uri)

        in_bucket = storage_client.bucket(in_bucket_name)
        in_blob = in_bucket.blob(in_blob_name)

        with tempfile.TemporaryDirectory() as tmpdir:
            local_fastq = os.path.join(tmpdir, os.path.basename(in_blob_name))
            in_blob.download_to_filename(local_fastq)

            metrics = compute_fastq_metrics(local_fastq)

            result = {
                **metrics,
                "input_uri": input_uri,
                "ran_at_utc": datetime.now(timezone.utc).isoformat(),
                "run_id": run_id,
            }

            output_name = f"qc-results/run_{run_id}.json"
            out_bucket = storage_client.bucket(results_bucket)
            out_blob = out_bucket.blob(output_name)
            out_blob.upload_from_string(
                json.dumps(result, indent=2),
                content_type="application/json",
            )

            output_uri = f"gs://{results_bucket}/{output_name}"

            mark_run_complete(
                conn=conn,
                run_id=run_id,
                read_count=result["reads"],
                total_bases=result["total_bases"],
                avg_read_length=result["avg_read_length"],
                gc_percent=result["gc_percent"],
                output_uri=output_uri,
            )

            print(json.dumps(result, indent=2))

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
