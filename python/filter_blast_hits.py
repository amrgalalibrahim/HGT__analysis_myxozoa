#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_blast_hits.py
====================
Filter tBLASTn results to retain only robust fungal-to-myxozoan candidates
for the HGT analysis.  Default thresholds match the Materials and Methods:

    identity    >= 60  (%)
    alignment   >= 100 (aa)
    E-value     <= 1e-30

Usage:
    python filter_blast_hits.py --blast_dir <dir> --out <file> \\
        --identity 60 --min_len 100 --max_evalue 1e-30

Inputs:  BLAST outfmt 6 files (*.tblastn) in --blast_dir
Outputs: single TSV with header

Author: Ibrahim A. G. A. <amrgalalibrahim@gmail.com>
"""
import argparse
import csv
import os
from pathlib import Path


FIELDS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--identity", type=float, default=60.0)
    ap.add_argument("--min_len", type=int, default=100)
    ap.add_argument("--max_evalue", type=float, default=1e-30)
    args = ap.parse_args()

    rows = []
    blast_dir = Path(args.blast_dir)
    for tblastn in sorted(blast_dir.glob("*.tblastn")):
        with open(tblastn) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                rec = dict(zip(FIELDS, parts[:12]))
                try:
                    pident = float(rec["pident"])
                    alen = int(rec["length"])
                    evalue = float(rec["evalue"])
                except ValueError:
                    continue
                if (pident >= args.identity and
                        alen >= args.min_len and
                        evalue <= args.max_evalue):
                    rec["source_file"] = tblastn.name
                    rows.append(rec)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS + ["source_file"],
                           dialect=csv.excel_tab)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {len(rows)} filtered hits to {out_path}")


if __name__ == "__main__":
    main()
