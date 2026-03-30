
   - demultiplexes reads by barcode
   - aligns reads with BWA
   - sorts and indexes BAMs with samtools
   - counts features with featureCounts
   - generates output artifacts
6. Outputs are written to the results bucket.
7. Run status and output metadata are updated in Cloud SQL.

## Current Inputs

### Ingest bucket
The pipeline scans a GCS folder for FASTQ files, for example:

`gs://<ingest-bucket>/data/`

### Marker file
A batch is triggered by uploading:

`gs://<ingest-bucket>/data/READY.txt`

### Reference assets
The analysis job expects reference files such as:

- `reference.fasta`
- `annotation.gtf`

These are typically stored in a separate reference bucket.

## Current Outputs

For each pipeline run, the job writes:

- **QC summary JSON**  
  Contains run-level statistics such as reads processed, assigned reads, GC percentage, and processed input files.

- **Count matrix CSV**  
  Gene-level counts by sample.

- **BED file**  
  Read-level genomic interval output derived from alignments.

Example output layout:

```text
gs://<results-bucket>/qc-results/run_<RUN_ID>.json
gs://<results-bucket>/count-matrices/count_matrix_run_<RUN_ID>.csv
gs://<results-bucket>/beds/reads_run_<RUN_ID>.bed


Barcode Demultiplexing

The current implementation uses fixed barcodes to split reads into samples before alignment.

Example mapping:

ACGTACGT → sample1
TGCATGCA → sample2
GATTACAG → sample3

These values are currently hard-coded in the analysis job and can be extended or externalized in future versions.

Tech Stack
Python 3.11
Google Cloud Storage
Google Cloud Functions Gen 2
Google Cloud Run Jobs
Google Cloud SQL
Docker
BWA
samtools
featureCounts
Repository Structure



Deployment Summary
Cloud Run Job

The Cloud Run Job runs the analysis container and expects environment variables such as:

DATA_DIR
REFERENCE_FASTA
ANNOTATION_GTF
RESULTS_BUCKET
RUN_ID
Trigger Function

The Cloud Function:

listens to the ingest bucket
ignores normal FASTQ uploads
launches the Cloud Run Job only when data/READY.txt is uploaded
Example Usage
Upload FASTQ files:
gsutil cp data/sample1.fastq gs://<ingest-bucket>/data/
gsutil cp data/sample2.fastq gs://<ingest-bucket>/data/
gsutil cp data/sample3.fastq gs://<ingest-bucket>/data/

Upload marker file:
echo ready > READY.txt
gsutil cp READY.txt gs://<ingest-bucket>/data/READY.txt

What This Project Demonstrates

This project is intended to show:

cloud-native pipeline orchestration
event-driven batch processing
managed container execution
run tracking with a relational database
integration of standard bioinformatics tools into a reproducible workflow
generation of structured output artifacts for downstream analysis
Current Status

This pipeline is a strong prototype and architecture demonstration. It uses standard alignment and counting tools, but current testing has focused primarily on synthetic datasets and workflow validation. Additional validation on public biological datasets, stricter filtering, and expanded QC are natural next steps.

Future Improvements

Planned or possible next steps include:

validation on real public datasets
stricter alignment filtering and QC thresholds
support for paired-end reads
external barcode/sample manifests
richer QC reporting
downstream visualization dashboards
workflow packaging with WDL or Nextflow
support for larger references and more realistic transcriptomic workflows
Why This Project Matters

Many pipelines work locally but fail when moved into a cloud production context. This project focuses on solving the engineering side of that transition: automated triggering, containerized execution, cloud-based tracking, and reproducible artifact generation.

It is designed as a practical example of how bioinformatics workflows can be deployed as managed, event-driven systems in the cloud.
