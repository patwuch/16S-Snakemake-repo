#!/usr/bin/env python3
"""Figures for the Diet x Radiation factorial diversity analysis (reza250909).

Ported from tisha251224's plot_diversity_factorial.py. Consumes what
run_diversity_factorial.R writes; only the PCoA ordinations are computed here,
from the same distance matrices, over exactly the samples the tests used.

Produces four figures, each as paired .svg + .png:
  1. alpha    - per-metric boxplots across the 2x2 design
  2. beta     - PCoA per distance metric, points with per-group ellipses
  3. permanova_r2 - PERMANOVA R2 per term, stacked per metric
  4. dispersion - betadisper distance-to-centroid per cell

Usage (from anywhere), inside the qiime2 container:
    python3 plot_diversity_factorial.py
"""
import argparse
import os
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
CORE = os.path.join(STUDY_DIR, "qiime2_artifacts", "core-metrics-results")
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
# Plots live in figures/, data and tables at the root. A hundred-odd plots at
# the root bury the handful of files a person actually opens.
FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)

DIET_ORDER = ["Control", "SCFAs"]
COND_ORDER = ["None", "IR"]
COND_MARKERS = {"None": "o", "IR": "^"}


def group_palette(meta_full):
    """Colour and order per Group, taken from the study's own metadata.

    Read from the UNSUBSETTED metadata on purpose: the same convention the main
    pipeline uses for GROUP_COLORS, so a group keeps its colour and position
    across every comparison, including ones it does not appear in.

    Group is the Diet x Radiation cell, so this also gives the cell palette.
    """
    for col in ("Group", "Diet", "Radiation", "PrimaryColor", "PrimaryOrder"):
        if col not in meta_full.columns:
            raise ValueError(f"metadata has no '{col}' column; cannot take colours from it")
    spec = meta_full[["Group", "PrimaryColor", "PrimaryOrder", "Diet", "Radiation"]].drop_duplicates()
    dup = spec.Group[spec.Group.duplicated()].unique()
    if len(dup):
        raise ValueError(f"Group(s) {list(dup)} carry more than one PrimaryColor/PrimaryOrder/cell")
    spec = spec.sort_values("PrimaryOrder")
    order = spec.Group.tolist()
    colors = dict(zip(spec.Group, spec.PrimaryColor))
    cell = {(d, r): g for g, d, r in zip(spec.Group, spec.Diet, spec.Radiation)}
    return order, colors, cell

ALPHA_PANELS = [("shannon", "any", "Shannon"), ("simpson", "any", "Simpson"),
                ("evenness", "any", "Pielou's Evenness"), ("chao1", "any", "Chao1"),
                ("faith_pd", "Greengenes2", "Faith's PD (Greengenes2)")]
BETA_PANELS = [("jaccard", "any", "jaccard-distance-matrix.qza", "Jaccard"),
               ("bray_curtis", "any", "bray-curtis-distance-matrix.qza", "Bray-Curtis"),
               ("unweighted_unifrac", "Greengenes2", "Greengenes2-unweighted-unifrac-distance-matrix.qza", "Unweighted UniFrac (Greengenes2)"),
               ("weighted_unifrac", "Greengenes2", "Greengenes2-weighted-unifrac-distance-matrix.qza", "Weighted UniFrac (Greengenes2)")]


def export_qza(qza, want):
    tmp = tempfile.mkdtemp()
    subprocess.run(["qiime", "tools", "export", "--input-path", qza, "--output-path", tmp],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    path = os.path.join(tmp, want)
    if not os.path.exists(path):
        raise FileNotFoundError(f"{want} not produced from {qza}")
    return path


def save(fig, stem):
    for ext, kw in (("svg", {}), ("png", {"dpi": 300})):
        path = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(path, format=ext, bbox_inches="tight", **kw)
        print(f"  wrote {path}")
    plt.close(fig)


def stat_label(stats, analysis, metric, db, term):
    row = stats[(stats.analysis == analysis) & (stats.metric == metric)
                & (stats.db == db) & (stats.term == term)]
    if row.empty or pd.isna(row.iloc[0]["p_raw"]):
        return f"{term}: n/a"
    p, q = row.iloc[0]["p_raw"], row.iloc[0]["p_BH"]
    # q is what the claim rests on; p is shown so a near-miss stays visible.
    mark = "*" if pd.notna(q) and q < 0.05 else ""
    return f"{term}: p={p:.3g}, q={q:.3g}{mark}"


def plot_alpha(alpha, stats, suffix, order, colors, cell):
    alpha = alpha.copy()
    alpha["Group"] = [cell[(d, r)] for d, r in zip(alpha.Diet, alpha.Radiation)]
    present = [g for g in order if g in set(alpha.Group)]
    palette = [colors[g] for g in present]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
    for ax in axes.ravel()[len(ALPHA_PANELS):]:
        ax.set_visible(False)
    for ax, (metric, db, label) in zip(axes.ravel(), ALPHA_PANELS):
        d = alpha[(alpha.metric == metric) & (alpha.db == db)]
        if d.empty:
            ax.set_visible(False)
            continue
        sns.boxplot(data=d, x="Group", y="value", order=present, palette=palette,
                    showfliers=False, linewidth=1.1, ax=ax)
        np.random.seed(42)   # stripplot's jitter draws from numpy's global RNG
        sns.stripplot(data=d, x="Group", y="value", order=present, size=5,
                      linewidth=0.6, edgecolor="black", color="white", ax=ax)
        # Headroom for the stats box, so it never sits on the highest points.
        lo, hi = d.value.min(), d.value.max()
        ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.35 * (hi - lo))
        ax.set_title(label, fontsize=11, weight="bold")
        ax.set_xlabel(""); ax.set_ylabel(label)
        ax.tick_params(axis="x", rotation=35)
        for lbl in ax.get_xticklabels():
            lbl.set_ha("right")
        sub = "\n".join(stat_label(stats, "alpha", metric, db, t)
                        for t in ("Diet", "Radiation", "Diet:Radiation"))
        ax.text(0.98, 0.98, sub, transform=ax.transAxes, fontsize=7.5,
                va="top", ha="right", family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))
    fig.suptitle("Alpha diversity, Diet x Radiation (n = 8 per cell)\n"
                 "two-way ANOVA on Diet * Radiation, Type II; q = BH within each term across the 5 metrics shown",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, f"diversity_alpha{suffix}")


def plot_beta(meta, suffix, order, colors, cell):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (metric, db, fname, label) in zip(axes.ravel(), BETA_PANELS):
        qza = os.path.join(CORE, fname)
        if not os.path.exists(qza):
            ax.set_visible(False)
            continue
        dm = pd.read_csv(export_qza(qza, "distance-matrix.tsv"), sep="\t", index_col=0)
        ids = [s for s in meta.index if s in dm.index]
        sub = dm.loc[ids, ids]
        ord_ = pcoa(DistanceMatrix(sub.values, ids=ids))
        pts = ord_.samples.iloc[:, :2]
        pts.index = ids
        ve = ord_.proportion_explained[:2] * 100
        d = meta.loc[ids].copy()
        d["PC1"], d["PC2"] = pts.iloc[:, 0].values, pts.iloc[:, 1].values

        for diet in DIET_ORDER:
            for cond in COND_ORDER:
                sel = d[(d.Diet == diet) & (d.Radiation == cond)]
                if sel.empty:
                    continue
                ax.scatter(sel.PC1, sel.PC2, c=colors[cell[(diet, cond)]],
                           marker=COND_MARKERS[cond], s=95, alpha=0.9,
                           edgecolor="black", linewidth=0.6, zorder=3)
                # Visual grouping only. A covariance ellipse from 8 points is not
                # a confidence region, and the report says so.
                if len(sel) >= 3:
                    cov = np.cov(sel.PC1, sel.PC2)
                    vals, vecs = np.linalg.eigh(cov)
                    theta = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
                    w, h = 2 * np.sqrt(np.maximum(vals[::-1], 0))
                    shade = colors[cell[(diet, cond)]]
                    ax.add_patch(Ellipse((sel.PC1.mean(), sel.PC2.mean()), w, h,
                                         angle=theta, facecolor=shade, alpha=0.10,
                                         edgecolor=shade, linestyle="--", linewidth=1.0,
                                         zorder=1))
        ax.set_title(label, fontsize=11, weight="bold")
        ax.set_xlabel(f"PC1 ({ve.iloc[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({ve.iloc[1]:.1f}%)")

    handles = [Line2D([], [], marker=COND_MARKERS[r], ls="", mfc=colors[cell[(d, r)]],
                      mec="black", ms=9, label=cell[(d, r)])
               for d in DIET_ORDER for r in COND_ORDER]
    handles.sort(key=lambda h: order.index(h.get_label()))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=10)
    # Deliberately bare: the report reads this plot on its own first, then brings
    # in PERMANOVA and PERMDISP, which run on the full distance matrix, not on these axes.
    fig.suptitle("Beta diversity, PCoA (first two axes)", fontsize=12)
    fig.tight_layout(rect=[0, 0.07, 1, 0.93])
    save(fig, f"diversity_beta_pcoa{suffix}")


# PERMANOVA terms are not groups, so they take a categorical palette of their own
# rather than the metadata group colours. Validated for colour-vision deficiency;
# the aqua slot is under 3:1 on white, which the direct labels cover.
TERM_COLORS = {"Diet": "#2a78d6", "Radiation": "#eb6834",
               "Diet:Radiation": "#1baf7a", "Residual": "#dcdcd7"}
TERM_LABELS = {"Diet": "Diet", "Radiation": "Radiation",
               "Diet:Radiation": "Interaction", "Residual": "Residual"}


def plot_permanova_r2(stats, suffix):
    """Stacked R2 per metric: how much of the variation in distances each term
    explains. The q-values cannot show effect size, and several sit at the
    permutation floor where they all look the same."""
    perm = stats[stats.analysis == "beta_permanova"]
    rows = []
    for metric, db, _f, label in BETA_PANELS:
        r = perm[(perm.metric == metric) & (perm.db == db)].set_index("term")["effect"]
        if r.empty:
            continue
        # Type II R2 plus the full-model residual need not sum to exactly 1 in
        # an unbalanced design; each segment is drawn at its own value.
        parts = {t: float(r[t]) for t in ("Diet", "Radiation", "Diet:Radiation", "Residual")}
        rows.append((label, parts))
    rows = rows[::-1]  # first metric on top

    fig, ax = plt.subplots(figsize=(10, 3.6))
    for y, (label, parts) in enumerate(rows):
        left = 0.0
        for term, val in parts.items():
            ax.barh(y, val * 100, left=left * 100, height=0.6, color=TERM_COLORS[term],
                    edgecolor="white", linewidth=2)
            if val >= 0.04:
                ax.text((left + val / 2) * 100, y, f"{val * 100:.0f}%", ha="center",
                        va="center", fontsize=9, color="#14181d")
            left += val
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.set_xlim(0, max(100, max(sum(p.values()) for _l, p in rows) * 100))
    ax.set_xlabel("share of variation in distances (PERMANOVA R², %)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=TERM_COLORS[t]) for t in TERM_COLORS]
    ax.legend(handles, [TERM_LABELS[t] for t in TERM_COLORS], ncol=4, frameon=False,
              loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=9)
    fig.tight_layout()
    save(fig, f"diversity_permanova_r2{suffix}")


def plot_dispersion(disp, stats, suffix, order, colors, cell):
    # betadisper labels cells "Diet/Radiation"; show them by their metadata
    # Group name, as every other figure does.
    disp = disp[disp.grouping == "Cell"].copy()
    disp["Group"] = [cell[tuple(c.split("/"))] for c in disp.group]
    present = [g for g in order if g in set(disp.Group)]
    palette = [colors[g] for g in present]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (metric, db, _f, label) in zip(axes.ravel(), BETA_PANELS):
        d = disp[(disp.metric == metric) & (disp.db == db)]
        if d.empty:
            ax.set_visible(False)
            continue
        sns.boxplot(data=d, x="Group", y="distance", order=present, palette=palette,
                    showfliers=False, linewidth=1.1, ax=ax)
        np.random.seed(42)
        sns.stripplot(data=d, x="Group", y="distance", order=present, size=5,
                      color="black", alpha=0.65, ax=ax)
        ax.set_title(label, fontsize=11, weight="bold")
        ax.set_xlabel(""); ax.set_ylabel("distance to centroid")
        ax.tick_params(axis="x", rotation=35)
        for lbl in ax.get_xticklabels():
            lbl.set_ha("right")
        ax.text(0.02, 0.97, stat_label(stats, "beta_dispersion", metric, db, "Cell"),
                transform=ax.transAxes, fontsize=8, va="top", ha="left",
                family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.85))
    fig.suptitle("Within-cell dispersion (PERMDISP distance to centroid)\n"
                 "spread, not location: tests whether a PERMANOVA result could be a "
                 "variance difference rather than a centroid shift", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save(fig, f"diversity_dispersion{suffix}")


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    suffix = ""

    # Read labels literally, so no level name is ever mistaken for a missing value.
    literal = dict(keep_default_na=False, na_values=[""])
    stats = pd.read_csv(os.path.join(TAB_DIR, f"diversity_factorial{suffix}.tsv"), sep="\t")
    alpha = pd.read_csv(os.path.join(TAB_DIR, f"diversity_alpha_values{suffix}.tsv"), sep="\t", **literal)
    disp = pd.read_csv(os.path.join(TAB_DIR, f"diversity_dispersion_distances{suffix}.tsv"), sep="\t", **literal)

    # Full metadata for the palette, subset only for the data being drawn.
    meta_full = pd.read_csv(os.path.join(OUT_DIR, "metadata_factorial.tsv"), sep="\t", index_col=0, **literal)
    order, colors, cell = group_palette(meta_full)
    meta = meta_full
    for col, expected in (("Diet", DIET_ORDER), ("Radiation", COND_ORDER)):
        missing = set(expected) - set(meta[col].unique())
        if missing:
            raise ValueError(f"metadata {col} is missing levels {missing}; got {sorted(meta[col].unique())}")

    print(f"[plot] {len(meta)} samples, {stats.metric.nunique()} metrics")
    plot_alpha(alpha, stats, suffix, order, colors, cell)
    plot_beta(meta, suffix, order, colors, cell)
    plot_permanova_r2(stats, suffix)
    plot_dispersion(disp, stats, suffix, order, colors, cell)


if __name__ == "__main__":
    main()
