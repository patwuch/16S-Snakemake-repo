#!/usr/bin/env python3
"""picrust2_nsti for reza250909's factorial report: weighted NSTI per sample, by
group (prediction quality), from run_picrust2.sh's output.

Usage, inside the qiime2 container:
    python3 plot_picrust2_factorial.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plot_diversity_factorial import group_palette, save

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
PIC = os.path.join(OUT_DIR, "picrust2", "pipeline")


def main():
    meta = pd.read_csv(os.path.join(OUT_DIR, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False)
    order, colors, _ = group_palette(meta)

    wn = pd.read_csv(os.path.join(PIC, "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0)
    d = meta.join(wn, how="inner")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(data=d, x="Group", y="weighted_NSTI", order=order, palette=[colors[g] for g in order],
                showfliers=False, linewidth=1.1, ax=ax)
    np.random.seed(42)   # stripplot's jitter draws from numpy's global RNG
    sns.stripplot(data=d, x="Group", y="weighted_NSTI", order=order, size=5, linewidth=0.6,
                  edgecolor="black", color="white", ax=ax)
    ax.set_ylim(0, max(0.2, d.weighted_NSTI.max() * 1.15))
    ax.set_xlabel(""); ax.set_ylabel("weighted NSTI")
    ax.set_title("Distance of each mouse's ASVs to their nearest reference genome", fontsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    save(fig, "picrust2_nsti")


if __name__ == "__main__":
    main()
