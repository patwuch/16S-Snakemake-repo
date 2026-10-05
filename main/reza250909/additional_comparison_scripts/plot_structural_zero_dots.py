#!/usr/bin/env python3
"""Dot matrix of the structural-zero taxa: how many mice carry each, and how much they carry.

The group-level companion to plot_untestable_heatmap.py, with the same rows in the same order.
Columns are the Diet x Radiation groups. A dot's area is the number of the group's mice that
carry the taxon; its colour is the carriers' median share of reads (log scale, from the share
one read is of the shallowest mouse). A group with no carrier has no dot, only a dashed ring:
that group is where the taxon is a structural zero.

Shares are of each mouse's unrarefied reads in study-seqs.biom, collapsed to species level as
run_ancombc2_factorial.R does, the same numbers as structural_zeros.py.

Usage, inside the qiime2 container:
    python3 plot_structural_zero_dots.py [--out DIR]     (default: the report's figures/)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from plot_interaction_forest import RANK_WORD, axis_labels
from plot_untestable_heatmap import MUTED_INK, RAMP, RULE_INK, TEXT_INK
from structural_zeros import META, OUT_DIR, STEM, TAB_DIR, abundance, one_read_pct, pct_str

MAX_AREA = 330   # points² for a dot carried by every mouse of a group
RING_AREA = 60   # the dashed ring of an empty group: small, so it does not read as a count


def rows(un, groups, cell_cols):
    """The heatmap's row order: blocks by the set of groups missing the taxon, one missing group
    before two, larger blocks first, ties in group order; most widespread first within a block."""
    un = un.copy()
    un["absent"] = [tuple(g for g in groups if r[cell_cols[g]] == 0) for _, r in un.iterrows()]
    size = un.absent.value_counts()
    un = un.assign(k1=[len(p) for p in un.absent], k2=[-size[p] for p in un.absent],
                   k3=[tuple(groups.index(g) for g in p) for p in un.absent], k4=-un.n_samples)
    return un.sort_values(["k1", "k2", "k3", "k4"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "figures"))
    out_dir = ap.parse_args().out

    meta, pct = abundance()
    spec = pd.read_csv(META, sep="\t", index_col=0, keep_default_na=False)[["Group", "PrimaryColor", "PrimaryOrder"]]
    spec = spec.drop_duplicates().sort_values("PrimaryOrder")
    groups, colors = spec.Group.tolist(), dict(zip(spec.Group, spec.PrimaryColor))
    un = pd.read_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t")
    cell_cols = {g: next(c for c in un.columns if c.startswith(f"{g} (n=")) for g in groups}
    un = rows(un, groups, cell_cols)

    n_mice = {g: int((meta.Group == g).sum()) for g in groups}
    carriers = np.zeros((len(un), len(groups)), dtype=int)
    median = np.full((len(un), len(groups)), np.nan)
    for j, g in enumerate(groups):
        v = pct.loc[un.taxon, meta.index[meta.Group == g]]
        carriers[:, j] = (v > 0).sum(axis=1).values
        median[:, j] = [np.median(r[r > 0]) if (r > 0).any() else np.nan for r in v.values]
        if not np.array_equal(carriers[:, j], un[cell_cols[g]].values):
            raise SystemExit(f"carriers in {g} disagree with the untestable table")

    floor = one_read_pct()
    cmap = LinearSegmentedColormap.from_list("seq", RAMP)
    norm = LogNorm(vmin=floor, vmax=np.nanmax(median))

    nrow, ncol = len(un), len(groups)
    fig, ax = plt.subplots(figsize=(7.6, 0.27 * nrow + 2.4))
    xs = np.arange(ncol) * 1.0
    for i in range(nrow):
        for j in range(ncol):
            k = carriers[i, j]
            if k:
                ax.scatter(xs[j], i, s=MAX_AREA * k / n_mice[groups[j]], color=cmap(norm(median[i, j])),
                           edgecolor="white", linewidth=0.4, zorder=3)
            else:
                ax.scatter(xs[j], i, s=RING_AREA, facecolor="none", edgecolor=MUTED_INK,
                           linewidth=0.7, linestyle=(0, (2, 1.6)), zorder=2)

    # Group band and names above the columns.
    for j, g in enumerate(groups):
        ax.add_patch(Rectangle((xs[j] - 0.42, -1.25), 0.84, 0.32, color=colors[g], lw=0, clip_on=False))
        ax.text(xs[j], -1.45, f"{g}\n(n = {n_mice[g]})", ha="center", va="bottom", fontsize=8.5, color=TEXT_INK)

    labels = [lab if "unresolved" in lab or not t.rstrip().endswith("__") else
              f"{lab} · {RANK_WORD.get(t.rstrip().split(';')[-1].strip()[:1], 'rank')} unresolved"
              for lab, t in zip(axis_labels(list(un.taxon)), un.taxon)]
    ax.set_yticks(np.arange(nrow), labels, fontsize=7.5, color=TEXT_INK)
    for tick in ax.get_yticklabels():
        tick.set_fontstyle("italic")

    # Blocks by absence pattern, as in the heatmap, labelled on the right.
    right = xs[-1] + 0.62
    start = 0
    for i in range(1, nrow + 1):
        if i == nrow or un.absent[i] != un.absent[start]:
            if start:
                ax.axhline(start - 0.5, color=MUTED_INK, lw=0.6, clip_on=False)
            gs = list(un.absent[start])
            names = gs[0] if len(gs) == 1 else ", ".join(gs[:-1]) + " and " + gs[-1]
            ax.text(right, (start + i - 1) / 2, "not in " + names, va="center", ha="left", fontsize=7.5, color=MUTED_INK)
            start = i
    for i in range(nrow):
        ax.axhline(i, color=RULE_INK, lw=0.4, zorder=1, xmin=0.02, xmax=0.98)

    ax.set_xlim(xs[0] - 0.55, xs[-1] + 0.55)
    ax.set_ylim(nrow - 0.5, -1.55)
    ax.set_xticks([])
    ax.tick_params(axis="y", length=0, pad=6)
    for s in ax.spines.values():
        s.set_visible(False)

    # Legends: dot size (carriers), the empty ring, colour (median share of carriers).
    ref = n_mice[groups[0]]
    handles = [Line2D([], [], ls="", marker="o", markersize=np.sqrt(MAX_AREA * k / ref), markerfacecolor=RAMP[2],
                      markeredgecolor="white", label=f"{k} of {ref}") for k in (1, 2, 4, 8)]
    handles.append(Line2D([], [], ls="", marker="o", markersize=np.sqrt(RING_AREA), markerfacecolor="none",
                          markeredgecolor=MUTED_INK, markeredgewidth=0.7, label="none: structural zero"))
    leg = fig.legend(handles=handles, title="mice carrying it", loc="lower left", bbox_to_anchor=(0.02, 0.005),
                     ncol=5, frameon=False, fontsize=7.5, title_fontsize=8, handletextpad=0.9, columnspacing=1.3,
                     borderpad=0.2)
    leg.get_title().set_color(TEXT_INK)
    leg._legend_box.sep = 9   # room between the title and the largest dot
    cax = fig.add_axes([0.72, 0.045, 0.24, 0.011])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation="horizontal")
    # No tick close to the one-read end, where the labels would collide.
    ticks = [t for t in (0.01, 0.1, 1) if 2.5 * norm.vmin <= t <= norm.vmax]
    cb.set_ticks([norm.vmin] + ticks, labels=[f"{pct_str(norm.vmin)}%\none read"] + [f"{t:g}%" for t in ticks])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=7, colors=MUTED_INK)
    cb.outline.set_visible(False)
    cb.set_label("carriers' median share of reads (log)", fontsize=8, color=TEXT_INK, labelpad=3)
    cb.ax.xaxis.set_label_position("top")

    fig.suptitle(f"The {nrow} structural-zero taxa, by group", fontsize=11.5, color=TEXT_INK, y=0.995)
    fig.text(0.5, 0.968, "Greengenes2, species level · dot size: mice that carry the taxon · colour: how much of "
             "their reads it is (median) · dashed ring: absent from every mouse of the group",
             ha="center", fontsize=7.5, color=MUTED_INK)
    fig.subplots_adjust(left=0.40, right=0.80, top=0.935, bottom=0.10)
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "svg"):
        out = os.path.join(out_dir, f"{STEM}_structural_zero_dots.{ext}")
        fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.15)
        print(f"  wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
