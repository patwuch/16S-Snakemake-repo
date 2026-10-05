#!/usr/bin/env python3
"""Weighted NSTI per sample, by group: how close each sample's ASVs sit to a reference genome
(run_picrust2.sh's KO_metagenome_out/weighted_nsti.tsv.gz). Writes figures/picrust2_nsti.

Ported from part 1 of reza250909/additional_comparison_scripts/plot_picrust2_factorial.py.

Usage, inside the qiime2 container:
    python3 plot_picrust2_nsti.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from plot_diversity_factorial import group_palette, save  # noqa: E402

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
PIC = os.path.join(OUT_DIR, "picrust2", "pipeline")


def main():
    meta = pd.read_csv(os.path.join(OUT_DIR, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False)
    order, colors, _cell = group_palette(meta)
    wn = pd.read_csv(os.path.join(PIC, "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0)
    d = meta.join(wn, how="inner")
    if len(d) != len(meta):
        raise SystemExit(f"weighted NSTI for {len(d)} of {len(meta)} samples")
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.boxplot(data=d, x="Group", y="weighted_NSTI", order=order, palette=[colors[g] for g in order],
                showfliers=False, linewidth=1.1, ax=ax)
    sns.stripplot(data=d, x="Group", y="weighted_NSTI", order=order, size=5, linewidth=0.6,
                  edgecolor="black", color="white", ax=ax)
    ax.set_ylim(0, max(0.2, d.weighted_NSTI.max() * 1.15))
    ax.set_xlabel(""); ax.set_ylabel("weighted NSTI")
    ax.set_title("Distance of each sample's ASVs to their nearest reference genome", fontsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    save(fig, "picrust2_nsti")


if __name__ == "__main__":
    main()
