import json
import os
from typing import Any, Dict

from google.cloud import run_v2

from db import create_pipeline_run, get_db_connection, upsert_sample


CONTAINER_NAME = os.getenv("JOB_CONTAINER_NAME", "worker")


def gcs_finalize_handler(event: Dict[str, Any], context: Any = None):
    bucket = event.get("bucket")
    file_name = event.get("name")

    if not bucket or not file_name:
        print("Missing bucket or file name in event payload.")
        return ("Bad Request", 400)

    # Only launch the batch job when the marker file is uploaded
    if file_name != "data/READY.txt":
        print(f"Skipping non-marker file: {file_name}")
        return ("Skipped", 200)

    input_uri = f"gs://{bucket}/{file_name}"
    pipeline_version = os.getenv("PIPELINE_VERSION", "v2-sql")

    conn = get_db_connection()
    try:
        sample_id = upsert_sample(
            conn,
            input_uri=input_uri,
            file_name=file_name,
            bucket_name=bucket,
        )
        run_id = create_pipeline_run(
            conn,
            sample_id=sample_id,
            status="received",
            pipeline_version=pipeline_version,
        )
    finally:
        conn.close()

    project = os.environ["PROJECT_ID"]
    region = os.environ["REGION"]
    job_name = os.environ["JOB_NAME"]

    client = run_v2.JobsClient()
    job_path = client.job_path(project, region, job_name)

    request = run_v2.RunJobRequest(
        name=job_path,
        overrides=run_v2.RunJobRequest.Overrides(
            container_overrides=[
                run_v2.RunJobRequest.Overrides.ContainerOverride(
                    name=CONTAINER_NAME,
                    env=[
                        run_v2.EnvVar(name="INPUT_URI", value=input_uri),
                        run_v2.EnvVar(name="RUN_ID", value=str(run_id)),
                    ],
                )
            ]
        ),
    )

    operation = client.run_job(request=request)

    print(
        json.dumps(
            {
                "message": "Job launched",
                "run_id": run_id,
                "input_uri": input_uri,
                "job_path": job_path,
                "container_name": CONTAINER_NAME,
                "operation": operation.operation.name,
            }
        )
    )

    return (
        json.dumps(
            {
                "message": "Job launched",
                "run_id": run_id,
                "input_uri": input_uri,
            }
        ),
        200,
    )  
