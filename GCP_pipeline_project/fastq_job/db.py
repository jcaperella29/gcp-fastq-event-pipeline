import os
from typing import Optional

import psycopg
from psycopg.rows import dict_row


def get_db_connection():
    database_url = os.environ["DATABASE_URL"]
    return psycopg.connect(database_url, row_factory=dict_row)


def mark_run_running(conn, run_id: int):
    with conn.cursor() as cur:
        cur.execute(
            """
            UPDATE pipeline_runs
            SET status = 'running',
                started_at = NOW()
            WHERE run_id = %s
            """,
            (run_id,),
        )
        conn.commit()


def mark_run_complete(
    conn,
    run_id: int,
    read_count: int,
    total_bases: int,
    avg_read_length: float,
    gc_percent: float,
    output_uri: str,
):
    with conn.cursor() as cur:
        cur.execute(
            """
            UPDATE pipeline_runs
            SET status = 'complete',
                read_count = %s,
                total_bases = %s,
                avg_read_length = %s,
                gc_percent = %s,
                output_uri = %s,
                finished_at = NOW()
            WHERE run_id = %s
            """,
            (
                read_count,
                total_bases,
                avg_read_length,
                gc_percent,
                output_uri,
                run_id,
            ),
        )
        conn.commit()


def mark_run_failed(conn, run_id: int, error_message: str):
    with conn.cursor() as cur:
        cur.execute(
            """
            UPDATE pipeline_runs
            SET status = 'failed',
                error_message = %s,
                finished_at = NOW()
            WHERE run_id = %s
            """,
            (error_message[:4000], run_id),
        )
        conn.commit()
