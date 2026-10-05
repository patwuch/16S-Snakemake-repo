#!/usr/bin/env python3
"""Presence dot map for the Day17/Day25 report: the factorial report's presence heatmap
(plot_aldex3_factorial.py, presence_heatmap) with the four arm x day groups as columns.

Rows are the 106 taxa the differential-abundance section analyses (present in >= 2 samples of an
arm), clustered on how many mice in each group carry them. Dot size = mice in the group carrying
the taxon (of 4), shade = the group's mean share of reads, a mouse without the taxon counting
as 0. No per-taxon marks: a group without the taxon shows as a missing dot. Day 17 columns are circles, Day 25 columns triangles. Drawn twice: all
taxa, and "_gplus" with only taxa named no further than genus, on one shared shade scale.

Usage, inside the qiime2 container:
    python3 plot_presence_day17_day25.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage  # noqa: E402

from plot_maaslin3_day17_day25 import CELL_COLOR, CELLS, MUTED_INK, TEXT_INK, genus_or_above, save  # noqa: E402
from structural_zeros_day17_day25 import label  # noqa: E402

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(THIS_DIR), "comparison-day17-day25", "tables")
STEM = "presence_day17day25_Greengenes2"


def presence_heatmap(genus_only, name):
    cnt = pd.read_csv(os.path.join(TAB, "collapsed_counts_day17day25_Greengenes2.tsv"), sep="\t", index_col=0)
    smp = pd.read_csv(os.path.join(TAB, "samples_day17day25.tsv"), sep="\t", index_col=0).loc[cnt.columns]
    cell = smp.Arm.map({"IFA": "IFA", "IFA_IR": "IFA+IR"}) + " " + smp.Timepoint
    depth = cnt.sum()  # the collapsed table keeps every taxon, so these are the full depths
    in_arm = pd.DataFrame({a: (cnt.loc[:, smp.Arm == a] > 0).sum(axis=1) >= 2 for a in ("IFA", "IFA_IR")})
    cnt = cnt[in_arm.any(axis=1)]
    pres = pd.DataFrame({g: (cnt.loc[:, cell == g] > 0).sum(axis=1) for g in CELLS})
    rel = cnt / depth * 100
    mean = pd.DataFrame({g: rel.loc[:, cell == g].mean(axis=1) for g in CELLS})
    # One shade scale for both views, set by all taxa, so a shade means the same in each.
    shown = mean.values[pres.values > 0]
    norm = LogNorm(vmin=shown.min(), vmax=shown.max())
    n_all = len(pres)
    if genus_only:
        pres = pres[[genus_or_above(t) for t in pres.index]]
    Z = linkage(pres.values, method="average", metric="euclidean")
    order = [pres.index[i] for i in leaves_list(Z)]
    pres, mean = pres.loc[order], mean.loc[order]
    n = len(order)
    fig = plt.figure(figsize=(8.5, 2.3 + 0.2 * n))
    h = fig.get_figheight()
    top, bot = 1 - 1.3 / h, 1.0 / h
    ax_d = fig.add_axes([0.01, bot, 0.08, top - bot])
    ax = fig.add_axes([0.10, bot, 0.24, top - bot])
    dendrogram(Z, ax=ax_d, orientation="left", no_labels=True, color_threshold=0, above_threshold_color=MUTED_INK)
    ax_d.invert_yaxis()
    ax_d.axis("off")
    # Sequential blue from step 200 up, log scale: the palest dot must stay visibly blue.
    cmap = LinearSegmentedColormap.from_list("seq", ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
    area = lambda k: 150 * k / 4  # dot area proportional to mice, 4 mice fills the cell
    for i in range(n):
        ax.axhline(i, color="#eeeeeb", lw=0.5, zorder=0)
    yy, xx = np.nonzero(pres.values)
    # Day 17 columns as circles, Day 25 columns as triangles (user, 2026-10-04); colour stays the shade.
    for day, mk in (("Day17", "o"), ("Day25", "^")):
        sel = np.array([CELLS[x].endswith(day) for x in xx], dtype=bool)
        ax.scatter(xx[sel], yy[sel], s=area(pres.values[yy[sel], xx[sel]]) * (1.25 if mk == "^" else 1),
                   c=mean.values[yy[sel], xx[sel]], cmap=cmap, norm=norm, marker=mk,
                   edgecolors="white", linewidths=0.4, zorder=2)
    ax.axvline(1.5, color=TEXT_INK, lw=1.2)
    ax.set_xlim(-0.5, len(CELLS) - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(range(len(CELLS)), CELLS, rotation=45, ha="left", fontsize=8.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", pad=16)
    for i, g in enumerate(CELLS):
        ax.add_patch(plt.Rectangle((i - 0.5, -1.25), 1, 0.6, color=CELL_COLOR[g.split()[0]], clip_on=False))
    # No per-taxon marks: a group without the taxon shows as a missing dot.
    ax.set_yticks(range(n), [label(t) for t in order], fontsize=8.5, color=TEXT_INK)
    ax.yaxis.tick_right()
    for tl in ax.get_yticklabels():
        tl.set_fontstyle("italic")
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)

    # Shade key, ends labelled with the actual extremes so neither reads as a cap.
    cax = fig.add_axes([0.10, bot - 0.40 / h, 0.30, 0.10 / h])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation="horizontal")
    lo, hi = norm.vmin, norm.vmax
    ticks = [t for t in (0.001, 0.01, 0.1, 1, 10)
             if lo < t < hi and abs(np.log10(t / lo)) > 0.6 and abs(np.log10(t / hi)) > 0.6]
    cb.set_ticks([lo] + ticks + [hi], labels=[f"{lo:.1g}%"] + [f"{t:g}%" for t in ticks] + [f"{hi:.2g}%"])
    cb.ax.minorticks_off()
    cb.set_label("group mean share of reads (log scale; absent mice count as 0)", fontsize=8.5, color=MUTED_INK)
    cb.ax.tick_params(labelsize=8, colors=MUTED_INK)
    cb.outline.set_visible(False)
    # Size key.
    kax = fig.add_axes([0.48, bot - 0.52 / h, 0.24, 0.30 / h])
    kax.scatter(range(4), [0] * 4, s=[area(k) for k in range(1, 5)], color=MUTED_INK, edgecolors="white", lw=0.4)
    for k in range(4):
        kax.text(k, -0.9, str(k + 1), ha="center", va="top", fontsize=8, color=MUTED_INK)
    kax.text(1.5, -2.1, "mice with the taxon, of 4 (no dot: none)", ha="center", va="top", fontsize=8.5, color=MUTED_INK)

    kax.set_xlim(-0.6, 3.6)
    kax.set_ylim(-2.5, 0.8)
    kax.axis("off")
    fig.text(0.60, bot - 0.78 / h, "● Day 17     ▲ Day 25", ha="center", va="top", fontsize=8.5, color=MUTED_INK)
    scope = "taxa named no further than genus" if genus_only else f"{n_all} taxa seen in at least 2 samples of an arm"
    fig.suptitle(f"Presence and abundance by group: the {n} {scope}" if genus_only
                 else f"Presence and abundance by group: the {scope}", x=0.01, y=1 - 0.2 / h, ha="left",
                 fontsize=11.5, color=TEXT_INK)
    save(fig, name)


if __name__ == "__main__":
    presence_heatmap(False, f"{STEM}_heatmap")
    presence_heatmap(True, f"{STEM}_gplus_heatmap")
