# Event-Driven FASTQ Processing Pipeline (GCP)

## 🚀 Overview

This project demonstrates a **cloud-native, event-driven bioinformatics pipeline** on Google Cloud Platform.

When a FASTQ file is uploaded to Cloud Storage, the system automatically:

1. Detects the upload event
2. Triggers a serverless function
3. Launches a containerized batch job
4. Computes basic FASTQ QC metrics
5. Writes results back to Cloud Storage

This mirrors real-world genomics data platform patterns.

---

## 🏗️ Architecture

```
Cloud Storage (ingest bucket)
        ↓ (object finalize event)
Eventarc Trigger
        ↓
Cloud Function (glue)
        ↓
Cloud Run Job (containerized compute)
        ↓
Cloud Storage (results bucket)
```

### Key properties

* Event-driven (no manual job submission)
* Serverless orchestration
* Containerized compute
* Cloud-native ETL pattern
* Easily extensible to real bioinformatics workflows

---

## 📦 Components

### 1. Cloud Function (`trigger_fn/`)

Responsible for:

* Receiving GCS finalize events
* Constructing the input URI
* Launching the Cloud Run Job

---

### 2. Cloud Run Job (`fastq_job/`)

Containerized batch job that:

* Downloads FASTQ from GCS
* Computes summary statistics:

  * read count
  * total bases
  * average read length
  * GC percentage
* Writes JSON summary to results bucket

---

### 3. Fake data generator

`make_fake_fastq.py` creates synthetic FASTQ data for testing.

---

## 🧪 Example output

```json
{
  "reads": 200,
  "total_bases": 16000,
  "avg_read_length": 80.0,
  "gc_percent": 49.8,
  "input_uri": "gs://.../sample.fastq",
  "ran_at_utc": "..."
}
```

---

## 🛠️ Tech Stack

* Google Cloud Functions (Gen2)
* Eventarc
* Cloud Run Jobs
* Cloud Storage
* Python 3.11
* Docker

---

## 🎯 Why this project matters

This project demonstrates patterns used in:

* Genomics data platforms
* Scientific ETL pipelines
* ML data ingestion systems
* Modern serverless architectures

It showcases:

* event-driven design
* cloud orchestration
* containerized batch processing
* production-style data flow

---

## 🔮 Possible Extensions

* Multi-stage workflows (Workflows DAG)
* BigQuery loading
* Parallel sample fan-out
* Retry/backoff logic
* Real FASTQ QC tools (FastQC, seqkit, etc.)
* Metadata tracking table
* Dashboard integration

---

## 👤 Author

JC — Bioinformatics / Scientific Software Engineering

---

## 📄 License

MIT
