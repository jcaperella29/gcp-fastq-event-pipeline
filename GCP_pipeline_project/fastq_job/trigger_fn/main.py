import os
from google.cloud import run_v2

REGION = os.environ.get("REGION", "us-central1")
JOB_NAME = os.environ["JOB_NAME"]           # e.g. fastq-summary-job
RESULTS_BUCKET = os.environ["RESULTS_BUCKET"]
PROJECT_ID = os.environ["PROJECT_ID"]

def gcs_to_job(event, context):
    # event has bucket + name
    bucket = event["bucket"]
    name = event["name"]
    print("=== FUNCTION INVOKED ===")
    print(event)
    # Ignore folder markers / empty
    if name.endswith("/"):
        return

    input_uri = f"gs://{bucket}/{name}"

    client = run_v2.JobsClient()
    parent = f"projects/{PROJECT_ID}/locations/{REGION}"
    job_path = f"{parent}/jobs/{JOB_NAME}"

    # Override env vars per execution (INPUT_URI changes per file)
    overrides = run_v2.RunJobRequest.Overrides(
        container_overrides=[
            run_v2.RunJobRequest.Overrides.ContainerOverride(
                env=[
                    run_v2.EnvVar(name="INPUT_URI", value=input_uri),
                    run_v2.EnvVar(name="RESULTS_BUCKET", value=RESULTS_BUCKET),
                ]
            )
        ]
    )

    req = run_v2.RunJobRequest(name=job_path, overrides=overrides)
    op = client.run_job(request=req)
    # Don't wait; just kick it off
    return f"Started job for {input_uri}"
