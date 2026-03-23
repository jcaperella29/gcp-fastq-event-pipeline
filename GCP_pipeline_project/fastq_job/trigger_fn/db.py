import os
from typing import Optional

import psycopg
from psycopg.rows import dict_row


def get_db_connection():
    """
    Basic Postgres connection.
    For local dev, use DATABASE_URL.
    Example:
      postgresql://user:password@host:5432/dbname
    """
    database_url = os.environ["DATABASE_URL"]
    return psycopg.connect(database_url, row_factory=dict_row)


def upsert_sample(conn, input_uri: str, file_name: str, bucket_name: str) -> int:
    with conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO samples (input_uri, file_name, bucket_name)
            VALUES (%s, %s, %s)
            ON CONFLICT (input_uri)
            DO UPDATE SET
                file_name = EXCLUDED.file_name,
                bucket_name = EXCLUDED.bucket_name
            RETURNING sample_id
            """,
            (input_uri, file_name, bucket_name),
        )
        row = cur.fetchone()
        conn.commit()
        return row["sample_id"]


def create_pipeline_run(
    conn,
    sample_id: int,
    status: str = "received",
    pipeline_version: Optional[str] = None,
) -> int:
    with conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO pipeline_runs (sample_id, status, pipeline_version)
            VALUES (%s, %s, %s)
            RETURNING run_id
            """,
            (sample_id, status, pipeline_version),
        )
        row = cur.fetchone()
        conn.commit()
        return row["run_id"]
