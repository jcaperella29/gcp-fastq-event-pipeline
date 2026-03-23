import json
import os
import subprocess
from typing import Any, Dict

from db import create_pipeline_run, get_db_connection, upsert_sample


def gcs_finalize_handler(event: Dict[str, Any], context: Any = None):
    """
    Triggered by GCS object finalize event.
    Expected event fields include:
      - bucket
      - name
    """

    bucket = event.get("bucket")
    file_name = event.get("name")

    if not bucket or not file_name:
        print("Missing bucket or file name in event payload.")
        return ("Bad Request", 400)

    # Optional: only process FASTQ-ish files
    if not (
        file_name.endswith(".fastq")
        or file_name.endswith(".fq")
        or file_name.endswith(".fastq.gz")
        or file_name.endswith(".fq.gz")
    ):
        print(f"Skipping non-FASTQ file: {file_name}")
        return ("Skipped", 200)

    input_uri = f"gs://{bucket}/{file_name}"
    pipeline_version = os.getenv("PIPELINE_VERSION", "v2-sql")

    conn = get_db_connection()
    try:
        sample_id = upsert_sample(conn, input_uri=input_uri, file_name=file_name, bucket_name=bucket)
        run_id = create_pipeline_run(
            conn,
            sample_id=sample_id,
            status="received",
            pipeline_version=pipeline_version,
        )
    finally:
        conn.close()

    project = os.environ["GCP_PROJECT"]
    region = os.environ["CLOUD_RUN_REGION"]
    job_name = os.environ["CLOUD_RUN_JOB_NAME"]

    # Launch Cloud Run Job and override env vars for this execution
    cmd = [
        "gcloud",
        "run",
        "jobs",
        "execute",
        job_name,
        "--region",
        region,
        "--project",
        project,
        "--update-env-vars",
        f"INPUT_URI={input_uri},RUN_ID={run_id}",
    ]

    print("Executing command:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    if result.returncode != 0:
        return (f"Failed to launch Cloud Run Job: {result.stderr}", 500)

    return (json.dumps({"message": "Job launched", "run_id": run_id, "input_uri": input_uri}), 200)
