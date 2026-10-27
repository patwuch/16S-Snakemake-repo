#!/usr/bin/env python3
"""Write tables/sample_reads.tsv: each mouse's total reads in study-seqs.biom
(unrarefied, every ASV), the library size run_picrust2_maaslin3.R reads.

Usage, inside the qiime2 container:
    python3 sample_reads.py
"""
import os

import biom
import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
META = os.path.join(OUT_DIR, "metadata_factorial.tsv")   # row order
OUT = os.path.join(OUT_DIR, "tables", "sample_reads.tsv")


def main():
    reads = biom.load_table(os.path.join(STUDY_DIR, "tables", "study-seqs.biom")).to_dataframe(dense=True).sum()
    order = pd.read_csv(META, sep="\t", index_col=0).index
    reads[order].astype(int).rename_axis("sample").rename("reads").to_csv(OUT, sep="\t")


if __name__ == "__main__":
    main()
