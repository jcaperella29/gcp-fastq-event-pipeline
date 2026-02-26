FROM python:3.11-slim

# system deps for bwa/samtools/subread(featureCounts)
RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa \
    samtools \
    subread \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY main.py /app/main.py

CMD ["python", "/app/main.py"]
