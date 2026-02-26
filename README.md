# Event-Driven FASTQ Processing Pipeline (GCP)

## 🚀 Overview

This project demonstrates a **cloud-native, event-driven bioinformatics pipeline** on Google Cloud Platform.

When a FASTQ file is uploaded to Cloud Storage, the system automatically:

1. Detects the upload event  
2. Triggers a serverless function  
3. Launches a containerized batch job  
4. Computes alignment and gene-level metrics  
5. Writes structured results to BigQuery  
6. Produces a pivoted gene × sample matrix  

This mirrors real-world genomics data platform patterns used in production environments.

---

## 🏗️ Architecture

Cloud Storage (ingest bucket)
↓ (object finalize event)
Eventarc Trigger
↓
Cloud Function (glue)
↓
Cloud Run Job (containerized compute)
↓
BigQuery (normalized tables)
↓
BigQuery SQL Pivot (gene × sample matrix)
Cloud Storage (ingest bucket)


### Key properties

* Event-driven (no manual job submission)
* Fully serverless orchestration
* Containerized bioinformatics compute
* Warehouse-ready normalized schema
* SQL-based matrix generation
* Easily extensible to real NGS workflows

---

## 📦 Components

### 1. Cloud Function (`trigger_fn/`)

Responsible for:

* Receiving GCS finalize events
* Constructing the FASTQ input URI
* Launching the Cloud Run Job with environment variables
* Acting as the event-driven control plane

---

### 2. Cloud Run Job (`fastq_job/`)

Containerized batch job that performs:

#### FASTQ processing

* Download FASTQ from GCS
* Align reads with **BWA**
* Process BAM with **samtools**
* Generate gene counts with **featureCounts**

#### QC + metrics

Computes:

* read count  
* total bases  
* average read length  
* GC percentage  

#### Warehouse loading

Writes structured results into BigQuery tables:

* `fastq_qc.qc_runs`
* `fastq_qc.contig_counts`
* `fastq_qc.gene_counts`

---

### 3. Fake data generator

`make_fake_fastq.py` creates synthetic FASTQ data for testing and pipeline validation.

---

## 🧪 Example BigQuery Pivot

After multiple samples are processed, the pipeline supports warehouse-native matrix generation:

```sql
SELECT *
FROM (
  SELECT sample_name, gene_id, count
  FROM `fastq-data-pipeline.fastq_qc.gene_counts`
  WHERE sample_name IS NOT NULL
)
PIVOT (
  SUM(count) FOR sample_name IN (
    "sample_align1.fastq" AS sample_align1,
    "sample_align2.fastq" AS sample_align2,
    "sample_align3.fastq" AS sample_align3
  )
)
ORDER BY gene_id;

🛠️ Tech Stack

Cloud & Orchestration

Google Cloud Functions (Gen2)

Eventarc

Cloud Run Jobs

Cloud Storage

BigQuery

Bioinformatics

BWA

samtools

featureCounts

Runtime

Python 3.11

Docker

🎯 Why this project matters

This project demonstrates patterns used in:

Genomics data platforms

Scientific ETL pipelines

ML feature pipelines

Production bioinformatics backends

Modern serverless architectures

It showcases:

event-driven design

cloud-native orchestration

containerized bioinformatics workloads

warehouse-ready data modeling

SQL-based analytical pivoting

end-to-end automated sample processing

🔮 Possible Extensions

Multi-sample fan-out

Workflows / Cromwell DAG integration

Partitioned BigQuery tables

Looker / dashboard layer

Retry & idempotency hardening

FastQC / seqkit integration

Metadata lineage tracking

Multi-region scaling

👤 Author

JC
Bioinformatics • Scientific Software Engineering • Cloud Genomics
