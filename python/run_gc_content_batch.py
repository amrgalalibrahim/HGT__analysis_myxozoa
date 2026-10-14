#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
	Horizontal Gene Transfer from Fungi to Myxozoa: Genomic Evidence and Evolutionary Implications		
		Amr Galal Abdelraheem Ibrahim; Edson Aparecido Adriano	

	Batch runner for GC content analysis across multiple genomes/candidate sets.

Input manifest TSV:
- sample_id
- genome_fasta
- candidates_tsv

Outputs:
- one gc_content_results.tsv per sample in --outdir
- a merged summary file: merged_gc_content.tsv

Python: 3.8+

Creator: Ibrahim A. G. A. (amrgalalibrahim@gmail.com)
"""

import argparse
import csv
import os
import subprocess
from typing import List

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GC_SCRIPT = os.path.join(SCRIPT_DIR, "gc_content_analysis.py")


def main():
    ap = argparse.ArgumentParser(description="Batch GC content analysis")
    ap.add_argument("--manifest", required=True, help="TSV with sample_id, genome_fasta, candidates_tsv")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--flank", type=int, default=10000, help="Flank size (bp)")
    ap.add_argument("--flag_delta", type=float, default=5.0, help="Flag threshold in pct points")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    merged_path = os.path.join(args.outdir, "merged_gc_content.tsv")
    merged_rows: List[dict] = []

    with open(args.manifest, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, dialect=csv.excel_tab)
        req = {"sample_id", "genome_fasta", "candidates_tsv"}
        if reader.fieldnames is None or not req.issubset(set(reader.fieldnames)):
            raise ValueError("manifest must include header: sample_id, genome_fasta, candidates_tsv")

        for row in reader:
            sid = row["sample_id"].strip()
            genome = row["genome_fasta"].strip()
            cand = row["candidates_tsv"].strip()
            out = os.path.join(args.outdir, f"{sid}_gc_content.tsv")

            cmd = [
                "python3", GC_SCRIPT,
                "--genome", genome,
                "--candidates", cand,
                "--flank", str(args.flank),
                "--flag_delta", str(args.flag_delta),
                "--out", out
            ]
            subprocess.check_call(cmd)

            # merge
            with open(out, "r", encoding="utf-8") as oh:
                r2 = csv.DictReader(oh, dialect=csv.excel_tab)
                for r in r2:
                    r["sample_id"] = sid
                    merged_rows.append(r)

    # write merged
    if merged_rows:
        fields = ["sample_id"] + [k for k in merged_rows[0].keys() if k != "sample_id"]
        with open(merged_path, "w", encoding="utf-8", newline="") as mh:
            w = csv.DictWriter(mh, fieldnames=fields, dialect=csv.excel_tab)
            w.writeheader()
            for r in merged_rows:
                w.writerow(r)


if __name__ == "__main__":
    main()
