#!/usr/bin/env python3
"""Figures for the PICRUSt2 section of reza250909's factorial report.

Reads run_picrust2.sh's output and run_picrust2_factorial.R's tables:
  1. picrust2_nsti            - weighted NSTI per sample, by group (prediction quality)
  2. picrust2_pathway_pcoa    - PCoA of Bray-Curtis distances between pathway profiles
  3. picrust2_pathway_effects - SCFAs+IR vs IR among pathways predicted in >= 16 of 32 mice:
                                every one with raw p < 0.05, log2FC with 95% CI

Usage, inside the qiime2 container:
    python3 plot_picrust2_factorial.py
"""
import html
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform

from plot_diversity_factorial import COND_MARKERS, DIET_ORDER, COND_ORDER, group_palette, save

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
TAB_DIR = os.path.join(OUT_DIR, "tables")
PIC = os.path.join(OUT_DIR, "picrust2", "pipeline")
INK, MUTED = "#14181d", "#6b7280"


def main():
    meta = pd.read_csv(os.path.join(OUT_DIR, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False)
    order, colors, cell = group_palette(meta)

    # 1. weighted NSTI per sample
    wn = pd.read_csv(os.path.join(PIC, "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0)
    d = meta.join(wn, how="inner")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(data=d, x="Group", y="weighted_NSTI", order=order, palette=[colors[g] for g in order],
                showfliers=False, linewidth=1.1, ax=ax)
    sns.stripplot(data=d, x="Group", y="weighted_NSTI", order=order, size=5, linewidth=0.6,
                  edgecolor="black", color="white", ax=ax)
    ax.set_ylim(0, max(0.2, d.weighted_NSTI.max() * 1.15))
    ax.set_xlabel(""); ax.set_ylabel("weighted NSTI")
    ax.set_title("Distance of each mouse's ASVs to their nearest reference genome", fontsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    save(fig, "picrust2_nsti")

    # 2. PCoA of pathway profiles
    ab = pd.read_csv(os.path.join(PIC, "pathways_out", "path_abun_unstrat.tsv.gz"), sep="\t", index_col=0)[meta.index]
    rel = ab / ab.sum()
    dm = squareform(pdist(rel.T.values, metric="braycurtis"))
    o = pcoa(DistanceMatrix(dm, ids=list(meta.index)))
    pts = o.samples.iloc[:, :2].values
    ve = o.proportion_explained[:2] * 100
    fig, ax = plt.subplots(figsize=(6, 5))
    for diet in DIET_ORDER:
        for cond in COND_ORDER:
            sel = ((meta.Diet == diet) & (meta.Radiation == cond)).values
            ax.scatter(pts[sel, 0], pts[sel, 1], c=colors[cell[(diet, cond)]], marker=COND_MARKERS[cond],
                       s=80, edgecolor="black", linewidth=0.6, zorder=3)
    ax.set_xlabel(f"PC1 ({ve.iloc[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ve.iloc[1]:.1f}%)")
    handles = [Line2D([], [], marker=COND_MARKERS[r], ls="", mfc=colors[cell[(dd, r)]], mec="black", ms=8,
                      label=cell[(dd, r)]) for dd in DIET_ORDER for r in COND_ORDER]
    handles.sort(key=lambda h: order.index(h.get_label()))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, frameon=False, fontsize=9)
    ax.set_title("Predicted MetaCyc pathway profiles, PCoA of Bray-Curtis distances", fontsize=10)
    sns.despine(ax=ax)
    fig.tight_layout()
    save(fig, "picrust2_pathway_pcoa")

    # 3. SCFAs+IR vs IR among well-predicted pathways (>= 16 of 32 mice), raw p < 0.05
    sw = pd.read_csv(os.path.join(TAB_DIR, "PICRUSt2_pathways_ancombc2_prevalence_sweep_long.tsv"), sep="\t")
    x = sw[(sw.min_presence == 16) & (sw.contrast == "SCFAs+IR vs IR") & (sw.p < 0.05)]
    n_tested = int(((sw.min_presence == 16) & (sw.contrast == "SCFAs+IR vs IR")).sum())
    x = x.sort_values("log2fc").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.5, 0.2 * len(x) + 1.6))
    y = np.arange(len(x))[::-1]
    col = np.where((x.log2fc > 0).values, colors["SCFAs+IR"], colors["IR"])
    ax.hlines(y, x.log2fc - 1.96 * x.se_log2, x.log2fc + 1.96 * x.se_log2, color=col, linewidth=1.2)
    ax.scatter(x.log2fc, y, facecolor="white", edgecolor=col, s=26, linewidth=1.3, zorder=3)
    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels([html.unescape(t) for t in x.description], fontsize=6.5, color=INK)
    ax.set_ylim(-1, len(x))
    ax.set_xlabel("log2 fold change, SCFAs+IR vs IR (95% CI)")
    ax.set_title(f"SCFAs+IR vs IR, {n_tested} pathways predicted in at least 16 of 32 mice:\n"
                 f"the {len(x)} with raw p < 0.05 (smallest q {x.q.min():.2f}; hollow = not a hit)", fontsize=9)
    handles = [Line2D([], [], marker="o", ls="", mfc="white", mec=colors[g], label=f"higher in {g}")
               for g, up in (("SCFAs+IR", True), ("IR", False)) if ((x.log2fc > 0) == up).any()]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=8)
    sns.despine(ax=ax)
    fig.tight_layout()
    save(fig, "picrust2_pathway_effects")


if __name__ == "__main__":
    main()
