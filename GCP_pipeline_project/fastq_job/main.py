import json
import os
import re
import sys
from datetime import datetime
from google.cloud import storage

def parse_gs_uri(gs_uri: str):
    m = re.match(r"^gs://([^/]+)/(.+)$", gs_uri)
    if not m:
        raise ValueError(f"Not a gs:// URI: {gs_uri}")
    return m.group(1), m.group(2)

def download_blob(gs_uri: str, local_path: str):
    bucket_name, blob_name = parse_gs_uri(gs_uri)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.download_to_filename(local_path)

def upload_text(bucket_name: str, blob_name: str, text: str, content_type="application/json"):
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.upload_from_string(text, content_type=content_type)

def fastq_summary(local_fastq: str):
    reads = 0
    total_len = 0
    gc = 0

    with open(local_fastq, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()

            if not (seq and plus and qual):
                break

            reads += 1
            total_len += len(seq)
            gc += sum(1 for c in seq if c in ("G", "C"))

    avg_len = (total_len / reads) if reads else 0.0
    gc_pct = (100.0 * gc / total_len) if total_len else 0.0

    return {
        "reads": reads,
        "total_bases": total_len,
        "avg_read_length": round(avg_len, 3),
        "gc_percent": round(gc_pct, 3),
    }

def main():
    # Prefer env vars (nice for orchestration)
    input_uri = os.environ.get("INPUT_URI") or (sys.argv[1] if len(sys.argv) > 1 else None)
    results_bucket = os.environ.get("RESULTS_BUCKET")

    if not input_uri or not results_bucket:
        raise SystemExit("Missing INPUT_URI and/or RESULTS_BUCKET")

    local_path = "/tmp/input.fastq"
    download_blob(input_uri, local_path)

    summary = fastq_summary(local_path)
    summary["input_uri"] = input_uri
    summary["ran_at_utc"] = datetime.utcnow().isoformat() + "Z"

    # Write results to gs://RESULTS_BUCKET/<original_name>.summary.json
    _, blob_name = parse_gs_uri(input_uri)
    out_name = blob_name.rsplit("/", 1)[-1] + ".summary.json"

    upload_text(results_bucket, out_name, json.dumps(summary, indent=2))
    print(f"Wrote gs://{results_bucket}/{out_name}")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()