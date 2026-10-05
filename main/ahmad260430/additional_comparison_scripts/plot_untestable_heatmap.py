#!/usr/bin/env python3
"""Heatmap of the taxa ANCOM-BC2 could not test because a whole group lacks them.

Rows are the features in tables/ANCOMBC2_factorial_Greengenes2_untestable_absent_from_a_cell.tsv
(in >= 2 samples of an arm but absent from every sample of at least one
Arm x Timepoint cell, written by build_untestable_day17_day25.py). Columns are
the 16 samples, grouped by cell, with each mouse in the same position at both
timepoints. A cell's colour is that feature's share of the sample's reads, on a
log scale; an empty cell means the feature was not seen in that sample.

Copied from reza250909/additional_comparison_scripts/plot_untestable_heatmap.py;
only the configuration block differs.

Counts are collapsed to the species-level Greengenes2 path exactly as
run_ancombc2_factorial.R collapses them (ancombc2_common.R, level 7), and
shares are of the sample's total reads in study-seqs.biom.

ahmad260430 has two analyses in one scripts folder, so this copy takes
--comparison to pick the configuration; the drawing code is the same.

Usage, inside the qiime2 container:
    python3 plot_untestable_heatmap.py [--comparison day17-day25|factorial]
"""
import argparse
import os

import biom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle

from plot_interaction_forest import RANK_WORD, axis_labels

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)

# ---- Study configuration: the only part that differs between the copies of
# this script in reza250909, tisha251224 and ahmad260430. ----
OUT_DIR = os.path.join(STUDY_DIR, "comparison-day17-day25")
META = os.path.join(OUT_DIR, "metadata_joint_cells.tsv")    # needs Group, PrimaryColor, PrimaryOrder
STEM = "DA_day17day25_Greengenes2"                         # <STEM>_untestable_absent_from_a_cell.tsv
UNIT = "samples"                                           # what a column is
RULE = "in ≥ 2 samples of an arm but absent from a whole arm × day cell" # why these rows were not tested
META_FILTER = None                                         # (column, value) to keep, or None

CONFIGS = {
    "day17-day25": {},
    "factorial": dict(
        OUT_DIR=os.path.join(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial"),
        META=os.path.join(STUDY_DIR, "metadata.tsv"),
        META_FILTER=("Timepoint", "Day17"),
        STEM="DA_factorial_Greengenes2_Day17_minpres6",
        UNIT="mice",
        RULE="seen in ≥ 6 mice but absent from a whole group"),
}
# ----

TAB_DIR = os.path.join(OUT_DIR, "tables")
FIG_DIR = os.path.join(OUT_DIR, "figures")
LEVEL = 7

# Sequential blue from step 200 up: the lightest present value must stay
# visibly blue, because white means absent.
RAMP = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
TEXT_INK, MUTED_INK, RULE_INK = "#2b2b28", "#6b6b66", "#dcdcd7"


def collapse(table, tax):
    """Sum ASV counts onto their first LEVEL taxonomy ranks, as collapse_to_level() does."""
    def path(t):
        if not isinstance(t, str):
            return None
        parts = [p.strip() for p in t.split(";") if p.strip()]
        return ";".join(parts[:LEVEL]) if parts else None
    paths = tax.reindex(table.index).map(path)
    return table[paths.notna()].groupby(paths[paths.notna()]).sum()


def configure(name):
    global OUT_DIR, META, META_FILTER, STEM, UNIT, RULE, TAB_DIR, FIG_DIR
    for k, v in CONFIGS[name].items():
        globals()[k] = v
    TAB_DIR = os.path.join(OUT_DIR, "tables")
    FIG_DIR = os.path.join(OUT_DIR, "figures")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison", choices=sorted(CONFIGS), default="day17-day25")
    configure(ap.parse_args().comparison)
    meta = pd.read_csv(META, sep="\t", index_col=0, keep_default_na=False, na_values=[""])
    if META_FILTER:
        meta = meta[meta[META_FILTER[0]] == META_FILTER[1]]
    spec = meta[["Group", "PrimaryColor", "PrimaryOrder"]].drop_duplicates().sort_values("PrimaryOrder")
    groups = spec.Group.tolist()
    colors = dict(zip(spec.Group, spec.PrimaryColor))

    counts = biom.load_table(os.path.join(STUDY_DIR, "tables", "study-seqs.biom")).to_dataframe(dense=True)
    counts = counts[meta.index]
    depth = counts.sum()
    tax = pd.read_csv(os.path.join(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                      sep="\t", index_col=0)["Taxon"]
    rel = collapse(counts, tax) / depth * 100

    un = pd.read_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t")
    missing = set(un.taxon) - set(rel.index)
    if missing:
        raise SystemExit(f"{len(missing)} untestable taxa not found after collapsing: {sorted(missing)[:3]}")
    cell_cols = {g: next(c for c in un.columns if c.startswith(f"{g} (n=")) for g in groups}
    # Check the recomputed presence against the table the model wrote.
    for g, c in cell_cols.items():
        mice = meta.index[meta.Group == g]
        got = (rel.loc[un.taxon, mice] > 0).sum(axis=1).values
        if not np.array_equal(got, un[c].values):
            raise SystemExit(f"presence in {g} disagrees with {c} in the untestable table")

    # Rows: one block per set of groups the taxon is missing from. Blocks run
    # from one missing group to two, larger blocks first, ties in group order;
    # within a block, most widespread taxon first.
    un["absent"] = [tuple(g for g in groups if r[cell_cols[g]] == 0) for _, r in un.iterrows()]
    pat_size = un.absent.value_counts()
    un = un.assign(k1=[len(p) for p in un.absent], k2=[-pat_size[p] for p in un.absent],
                   k3=[tuple(groups.index(g) for g in p) for p in un.absent], k4=-un.n_samples)
    un = un.sort_values(["k1", "k2", "k3", "k4"]).reset_index(drop=True)

    cols = [s for g in groups for s in meta.index[meta.Group == g]]
    data = rel.loc[un.taxon, cols].values
    shown = np.ma.masked_less_equal(data, 0)

    gap = 0.6  # blank columns between groups
    xpos, x = [], 0.0
    for g in groups:
        n = int((meta.Group == g).sum())
        xpos += [x + i for i in range(n)]
        x += n + gap

    nrow = len(un)
    fig, ax = plt.subplots(figsize=(11.5, 0.24 * nrow + 2.6))
    cmap = LinearSegmentedColormap.from_list("seq", RAMP)
    norm = LogNorm(vmin=max(shown.min(), 1e-3), vmax=shown.max())
    for j, xj in enumerate(xpos):
        for i in range(nrow):
            v = data[i, j]
            if v > 0:
                ax.add_patch(Rectangle((xj, i), 0.92, 0.84, color=cmap(norm(max(v, norm.vmin))), lw=0))
            else:
                ax.add_patch(Rectangle((xj, i), 0.92, 0.84, facecolor="white", edgecolor=RULE_INK, lw=0.5))

    # Group band and labels above the columns.
    x = 0.0
    for g in groups:
        n = int((meta.Group == g).sum())
        ax.add_patch(Rectangle((x, -1.05), n - 0.08, 0.35, color=colors[g], lw=0, clip_on=False))
        ax.text(x + (n - 0.08) / 2, -1.3, f"{g} (n = {n})", ha="center", va="bottom", fontsize=9, color=TEXT_INK)
        x += n + gap

    # Separators between absence patterns, labelled on the right.
    # A path that keeps an empty trailing rank (".;g__X;s__") is a separate
    # feature from one that stops at the genus (".;g__X"), so say which.
    labels = [lab if "unresolved" in lab or not t.rstrip().endswith("__") else
              f"{lab} · {RANK_WORD.get(t.rstrip().split(';')[-1].strip()[:1], 'rank')} unresolved"
              for lab, t in zip(axis_labels(list(un.taxon)), un.taxon)]
    ax.set_yticks(np.arange(nrow) + 0.42, labels, fontsize=7.5, color=TEXT_INK)
    for tick in ax.get_yticklabels():
        tick.set_fontstyle("italic")
    right = xpos[-1] + 1.3
    start = 0
    for i in range(1, nrow + 1):
        if i == nrow or un.absent[i] != un.absent[start]:
            if start:
                ax.axhline(start - 0.08, color=MUTED_INK, lw=0.6, xmin=0, xmax=1, clip_on=False)
            gs = list(un.absent[start])
            names = gs[0] if len(gs) == 1 else ", ".join(gs[:-1]) + " and " + gs[-1]
            ax.text(right, (start + i) / 2 - 0.08, "not in " + names,
                    va="center", ha="left", fontsize=7.5, color=MUTED_INK)
            start = i
    for i, n in enumerate(un.n_samples):
        ax.text(-0.5, i + 0.42, f"{n}", ha="right", va="center", fontsize=6.5, color=MUTED_INK)
    ax.text(-0.5, -0.45, UNIT, ha="right", va="center", fontsize=6.5, color=MUTED_INK)

    ax.set_xlim(-0.2, xpos[-1] + 1)
    ax.set_ylim(nrow, -1.6)
    ax.set_xticks([xj + 0.46 for xj in xpos], cols, rotation=90, fontsize=6.5, color=MUTED_INK)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", length=0, pad=22)
    for s in ax.spines.values():
        s.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.025, pad=0.06, aspect=45, shrink=0.55)
    one = UNIT.rstrip("s") if UNIT != "mice" else "mouse"
    cb.set_label(f"share of the {one}'s reads (%, log scale) · white = not seen in that {one}",
                 fontsize=8, color=TEXT_INK)
    # Half-decade ticks, and both ends labelled with the actual smallest and
    # largest share, so the bar's end is not read as a cap at the last round
    # tick. Round ticks within 0.2 decades of an end would collide and go.
    lo, hi = norm.vmin, norm.vmax
    ticks = [t for t in (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
             if abs(np.log10(t / lo)) > 0.5 and abs(np.log10(t / hi)) > 0.5 and lo < t < hi]
    labels = [f"{t:g}%" for t in ticks]
    cb.set_ticks([lo] + ticks + [hi], labels=[f"min {lo:.2g}%"] + labels + [f"max {hi:.3g}%"])
    cb.ax.minorticks_off()
    cb.ax.tick_params(labelsize=7, colors=MUTED_INK)
    cb.outline.set_visible(False)

    share = rel.loc[un.taxon].sum()
    # Per-column shares, for the report text and check_report.py.
    top = rel.loc[un.taxon].idxmax()
    pd.DataFrame({"sampleid": share.index, "Group": meta.loc[share.index, "Group"].values,
                  "untestable_pct": share.values, "top_taxon": top.values,
                  "top_taxon_pct": [rel.loc[t, m] for t, m in zip(top, share.index)]}
                 ).to_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_shares.tsv"),
                          sep="\t", index=False)
    fig.suptitle(f"The {nrow} structural-zero taxa: {RULE}",
                 fontsize=11.5, color=TEXT_INK, y=0.995)
    fig.text(0.5, 0.965, f"Greengenes2, species level · together {share.mean():.2f}% of an average {one}'s reads "
             f"(at most {share.max():.1f}%) · one column per {one}",
             ha="center", fontsize=8, color=MUTED_INK)
    fig.subplots_adjust(left=0.27, right=0.86, top=0.90, bottom=0.08)
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "svg"):
        out = os.path.join(FIG_DIR, f"{STEM}_untestable_heatmap.{ext}")
        fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.15)
        print(f"  wrote {os.path.basename(out)}")
    plt.close(fig)
    print(f"  {nrow} taxa · mean share {share.mean():.3f}% · max share {share.max():.3f}% "
          f"· max single cell {data.max():.3f}% · patterns {pat_size.to_dict()}")


if __name__ == "__main__":
    main()
