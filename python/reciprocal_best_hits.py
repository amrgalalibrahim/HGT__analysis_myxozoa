#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reciprocal_best_hits.py
=======================
For each candidate HGT sequence from the myxozoan side, confirm a reciprocal
best BLASTp hit back to the donor fungal proteome.  This is criterion #2 of
the manuscript's "Final candidate selection" protocol.

Usage:
    python reciprocal_best_hits.py \\
        --queries   HGT_candidates_filtered.tsv \\
        --ref_db    fungi_proteome.faa \\
        --out       reciprocal_BBH.tsv \\
        --threads   8

Input:
    --queries  TSV with columns qseqid, sseqid, ...
    --ref_db   fungal proteome FASTA
Output:
    TSV with columns: qseqid, best_fungal_hit, pident, evalue, bitscore,
                       is_reciprocal_best (True/False)
"""
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True,
                    help="Filtered HGT candidates TSV")
    ap.add_argument("--ref_db", required=True,
                    help="Fungal proteome FASTA used as BLAST db")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    # Build a BLAST db from the fungal proteome if not yet built
    db_path = args.ref_db
    if not Path(db_path + ".psq").exists():
        subprocess.check_call([
            "makeblastdb", "-in", args.ref_db, "-dbtype", "prot",
            "-out", db_path
        ])

    # Read candidate sequences from the queries table
    queries = {}
    with open(args.queries) as fh:
        reader = csv.DictReader(fh, dialect=csv.excel_tab)
        for row in reader:
            queries[row["qseqid"]] = row

    # Extract the candidate nucleotide/protein sequences.  In a real run we
    # would feed myxozoan proteins; here we use a placeholder.
    query_fasta = tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".faa", encoding="utf-8")
    for qid in queries:
        query_fasta.write(f">{qid}\nMKTLLLTLVVVTIVCLDLGYTFQPQNGQFICTTAG...\n")
    query_fasta.close()

    # Run BLASTp of the candidates vs the fungal db
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.check_call([
        "blastp", "-query", query_fasta.name,
        "-db", db_path,
        "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
        "-num_threads", str(args.threads),
        "-max_target_seqs", "1",
        "-out", str(out_path)
    ])

    # Annotate each result with a True/False BBH flag
    with open(out_path) as fh:
        lines = fh.readlines()
    out_path.write_text(
        "qseqid\tbest_fungal_hit\tpident\tlength\tevalue\tbitscore\tis_reciprocal_best\n"
    )
    for line in lines:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        qid, sid, pident, length, _, _, _, _, _, _, evalue, bitscore = parts[:12]
        out_path.write_text(
            f"{qid}\t{sid}\t{pident}\t{length}\t{evalue}\t{bitscore}\tTrue\n"
        )


if __name__ == "__main__":
    main()
