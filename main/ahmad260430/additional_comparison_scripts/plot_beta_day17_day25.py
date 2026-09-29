#!/usr/bin/env python3
"""Beta diversity figures for the Day17-vs-Day25 comparison, one arm at a time.

Consumes what run_diversity_day17_day25.py writes. The PCoAs are recomputed on
each arm's eight samples rather than subset from the stored 32-sample
ordination, whose axes were fitted with every other arm in the matrix.

Produces, each as .svg + .png in comparison-day17-day25/figures:
  diversity_day17day25_beta_pcoa        - 2 arms x 4 metrics, each mouse's two samples joined
  diversity_day17day25_permanova_r2     - Timepoint R2 per arm and metric
  diversity_day17day25_dispersion       - distance to timepoint median, one panel per metric

Usage, inside the qiime2 container:
    python3 plot_beta_day17_day25.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from qiime2 import Artifact
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
CORE = os.path.join(STUDY_DIR, "qiime2_artifacts", "core-metrics-results")
OUT_DIR = os.path.join(STUDY_DIR, "comparison-day17-day25")
FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")

ARMS = [("IFA", "IFA"), ("IFA+IR", "IFA_IR")]
DAY_ORDER = ["Day17", "Day25"]
DAY_COLORS = {"Day17": "#E69F00", "Day25": "#56B4E9"}  # same Okabe-Ito pair as run_comparison.py
PAIR_INK = "#9a9a97"
METRICS = [("jaccard", "db-independent", "jaccard-distance-matrix.qza", "Jaccard"),
           ("braycurtis", "db-independent", "bray-curtis-distance-matrix.qza", "Bray-Curtis"),
           ("unweighted_unifrac", "Greengenes2", "Greengenes2-unweighted-unifrac-distance-matrix.qza",
            "Unweighted UniFrac (GG2)"),
           ("weighted_unifrac", "Greengenes2", "Greengenes2-weighted-unifrac-distance-matrix.qza",
            "Weighted UniFrac (GG2)")]


def save(fig, stem):
    for ext, kw in (("svg", {}), ("png", {"dpi": 300})):
        path = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(path, format=ext, bbox_inches="tight", facecolor="white", **kw)
        print(f"  wrote {path}")
    plt.close(fig)


def join_mice(ax, frame, x, y):
    for _, pair in frame.groupby("MouseID"):
        pair = pair.set_index("Timepoint").reindex(DAY_ORDER)
        ax.plot(pair[x], pair[y], color=PAIR_INK, lw=0.9, zorder=1)


def day_handles():
    return [Line2D([], [], marker="o", ls="", mfc=DAY_COLORS[d], mec="black", ms=8, label=d)
            for d in DAY_ORDER] + [Line2D([], [], color=PAIR_INK, lw=0.9, label="same mouse")]


def plot_pcoa(metas):
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6))
    for row, (arm, slug) in enumerate(ARMS):
        md = metas[slug]
        for ax, (_m, _db, fname, label) in zip(axes[row], METRICS):
            dm = Artifact.load(os.path.join(CORE, fname)).view(DistanceMatrix).filter(md.index.tolist())
            ordn = pcoa(dm)
            d = md.loc[list(dm.ids)].copy()
            d["PC1"] = ordn.samples.iloc[:, 0].values
            d["PC2"] = ordn.samples.iloc[:, 1].values
            join_mice(ax, d, "PC1", "PC2")
            for day in DAY_ORDER:
                sel = d[d.Timepoint == day]
                ax.scatter(sel.PC1, sel.PC2, c=DAY_COLORS[day], s=80, edgecolor="black",
                           linewidth=0.6, zorder=3)
            ve = ordn.proportion_explained.iloc[:2] * 100
            ax.set_xlabel(f"PC1 ({ve.iloc[0]:.0f}%)", fontsize=9)
            ax.set_ylabel(f"PC2 ({ve.iloc[1]:.0f}%)", fontsize=9)
            ax.set_title(f"{arm} · {label}", fontsize=10.5, weight="bold")
            ax.tick_params(labelsize=7.5)
    fig.legend(handles=day_handles(), loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.suptitle("Beta diversity, PCoA fitted on each arm's eight samples (first two axes)", fontsize=12)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    save(fig, "diversity_day17day25_beta_pcoa")


def plot_r2(beta):
    rows = []
    for arm, _slug in ARMS:
        for metric, db, _f, label in METRICS:
            r = beta[(beta.arm == arm) & (beta.metric == metric) & (beta.db == db)]
            rows.append((f"{arm} · {label}", float(r.R2.iloc[0])))
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(10, 4.2))
    for y, (label, r2) in enumerate(rows):
        ax.barh(y, r2 * 100, height=0.6, color="#2a78d6", edgecolor="white", linewidth=2)
        ax.barh(y, (1 - r2) * 100, left=r2 * 100, height=0.6, color="#dcdcd7",
                edgecolor="white", linewidth=2)
        ax.text(r2 * 50, y, f"{r2 * 100:.0f}%", ha="center", va="center", fontsize=9, color="white")
    ax.axhline(3.5, color="#9a9a97", lw=0.8)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows], fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of variation in distances (PERMANOVA R², %)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#2a78d6", "#dcdcd7")],
              ["Timepoint", "Residual (between mice, within timepoint)"], ncol=2, frameon=False,
              loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=9)
    fig.tight_layout()
    save(fig, "diversity_day17day25_permanova_r2")


def plot_dispersion(disp, beta):
    # Laid out like the factorial report's dispersion figure: one panel per
    # metric, groups along x, box plus points, the test in a corner box. Each
    # mouse's two distances are joined within its arm.
    groups = [(arm, day) for arm, _slug in ARMS for day in DAY_ORDER]
    xpos = {g: i for i, g in enumerate(groups)}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (metric, db, _f, label) in zip(axes.ravel(), METRICS):
        d = disp[(disp.metric == metric) & (disp.db == db)].copy()
        d["x"] = [xpos[(a, t)] for a, t in zip(d.arm, d.Timepoint)]
        box = ax.boxplot([d.distance[d.x == i] for i in range(len(groups))], positions=range(len(groups)),
                         widths=0.6, showfliers=False, patch_artist=True,
                         medianprops=dict(color="black", lw=1.1),
                         boxprops=dict(lw=1.1), whiskerprops=dict(lw=1.1), capprops=dict(lw=1.1))
        for patch, (_arm, day) in zip(box["boxes"], groups):
            patch.set_facecolor(DAY_COLORS[day])
        for arm, _slug in ARMS:
            join_mice(ax, d[d.arm == arm], "x", "distance")
        ax.scatter(d.x, d.distance, s=22, color="black", alpha=0.65, zorder=3)
        ax.set_xticks(range(len(groups)), [f"{arm} {day}" for arm, day in groups],
                      rotation=35, ha="right")
        ax.axvline(1.5, color="0.8", lw=0.8)
        # Headroom so the stats box sits above the data.
        ax.set_ylim(0, d.distance.max() * 1.3)
        ax.set_title(label, fontsize=11, weight="bold")
        ax.set_ylabel("distance to timepoint median")
        lines = []
        for arm, _slug in ARMS:
            r = beta[(beta.arm == arm) & (beta.metric == metric) & (beta.db == db)].iloc[0]
            mark = "*" if r.permdisp_q < 0.05 else ""
            lines.append(f"{arm + ':':<8}p={r.permdisp_p:.3g}, q={r.permdisp_q:.3g}{mark}")
        ax.text(0.02, 0.97, "\n".join(lines), transform=ax.transAxes, fontsize=8, va="top",
                ha="left", family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.85))
    fig.legend(handles=day_handles(), loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.suptitle("Within-timepoint dispersion (PERMDISP distance to spatial median)\n"
                 "spread, not location: tests whether a PERMANOVA result could be a "
                 "variance difference rather than a shift", fontsize=12)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    save(fig, "diversity_day17day25_dispersion")


def main():
    metas = {}
    for _arm, slug in ARMS:
        md = pd.read_csv(os.path.join(TAB_DIR, f"metadata_subset_{slug}.tsv"), sep="\t", index_col=0)
        disp_ids = pd.read_csv(os.path.join(TAB_DIR, "diversity_day17day25_dispersion_distances.tsv"),
                               sep="\t").drop_duplicates("sample").set_index("sample").MouseID
        md["MouseID"] = disp_ids.reindex(md.index).values
        metas[slug] = md
    beta = pd.read_csv(os.path.join(TAB_DIR, "diversity_day17day25_beta_stats.tsv"), sep="\t")
    disp = pd.read_csv(os.path.join(TAB_DIR, "diversity_day17day25_dispersion_distances.tsv"), sep="\t")
    plot_pcoa(metas)
    plot_r2(beta)
    plot_dispersion(disp, beta)


if __name__ == "__main__":
    main()
