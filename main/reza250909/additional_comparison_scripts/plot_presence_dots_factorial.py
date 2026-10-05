#!/usr/bin/env python3
"""Clustered dot map of presence and abundance by group, as in ahmad260430's factorial report.

Rows are the 135 taxa the DA section covers: the 104 ANCOM-BC2 tests (tables/<STEM>.tsv) and the
31 structural zeros (<STEM>_untestable_absent_from_a_cell.tsv). Columns are the four groups, each
diet's unirradiated and irradiated groups side by side. Dot size is the number of the group's mice
that carry the taxon, and the shade is the group's mean share of reads, a mouse without the taxon
counting as 0; a group with no carrier has no dot. Rows are clustered on presence.

Drawn twice: every taxon, and "_gplus" with only the taxa named no further than genus, on the
same shade scale (a filter, not a refit). A sidecar .tsv per view gives the counts the report
quotes.

Usage, inside the qiime2 container:
    python3 plot_presence_dots_factorial.py [--out DIR]     (default: the report's figures/)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage  # noqa: E402

from structural_zeros import META, OUT_DIR, STEM, TAB_DIR, abundance  # noqa: E402

GROUPS = ["Control", "IR", "SCFAs", "SCFAs+IR"]   # each diet's radiation pair adjacent
TEXT_INK, MUTED_INK = "#2b2b28", "#6b6b66"
# Sequential blue from step 200 up, log scale: the palest dot must stay visibly blue.
RAMP = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def ranks(taxon):
    return [p.strip() for p in taxon.split(";") if p.strip() and not p.strip().endswith("__")]


def genus_or_above(taxon):
    return not ranks(taxon)[-1].startswith("s__")


def label(taxon):
    """Deepest rank assigned, worded as the heat trees word placeholders."""
    r = ranks(taxon)
    name = r[-1][3:]
    if r[-1].startswith("s__"):
        return name
    rank = {"g": "genus", "f": "family", "o": "order", "c": "class", "p": "phylum"}.get(r[-1][0], "higher rank")
    if taxon.rstrip().endswith("s__"):
        return f"{name} (no species call)" if rank == "genus" else f"{name} ({rank} only, no species call)"
    return f"{name} (classifier stopped at {rank})"


def draw(genus_only, pres, mean, norm, colors, n_mice, out_dir):
    if genus_only:
        keep = [genus_or_above(t) for t in pres.index]
        pres, mean = pres[keep], mean[keep]
    Z = linkage(pres.values, method="average", metric="euclidean")
    order = [pres.index[i] for i in leaves_list(Z)]
    pres, mean = pres.loc[order], mean.loc[order]
    n, k = len(order), n_mice
    fig = plt.figure(figsize=(8.2, 2.3 + 0.2 * n))
    h = fig.get_figheight()
    top, bot = 1 - 1.3 / h, 1.0 / h
    ax_d = fig.add_axes([0.01, bot, 0.08, top - bot])
    ax = fig.add_axes([0.10, bot, 0.26, top - bot])
    dendrogram(Z, ax=ax_d, orientation="left", no_labels=True, color_threshold=0, above_threshold_color=MUTED_INK)
    ax_d.invert_yaxis()
    ax_d.axis("off")
    cmap = LinearSegmentedColormap.from_list("seq", RAMP)
    area = lambda c: 150 * c / k   # dot area proportional to mice; every mouse fills the cell
    for i in range(n):
        ax.axhline(i, color="#eeeeeb", lw=0.5, zorder=0)
    yy, xx = np.nonzero(pres.values)
    ax.scatter(xx, yy, s=area(pres.values[yy, xx]), c=mean.values[yy, xx], cmap=cmap, norm=norm,
               edgecolors="white", linewidths=0.4, zorder=2)
    ax.axvline(1.5, color=TEXT_INK, lw=1.2)
    ax.set_xlim(-0.5, len(GROUPS) - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks(range(len(GROUPS)), GROUPS, rotation=45, ha="left", fontsize=8.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", pad=16)
    for i, g in enumerate(GROUPS):
        ax.add_patch(plt.Rectangle((i - 0.5, -1.25), 1, 0.6, color=colors[g], clip_on=False))
    # No per-taxon marks: a group without the taxon shows as a missing dot.
    ax.set_yticks(range(n), [label(t) for t in order], fontsize=8.5, color=TEXT_INK)
    ax.yaxis.tick_right()
    for tl in ax.get_yticklabels():
        tl.set_fontstyle("italic")
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)

    # Shade key, ends labelled with the actual extremes so neither reads as a cap.
    cax = fig.add_axes([0.10, bot - 0.40 / h, 0.26, 0.10 / h])
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
    steps = [1, 2, 4, 8] if k == 8 else list(range(1, k + 1))
    kax = fig.add_axes([0.44, bot - 0.52 / h, 0.24, 0.30 / h])
    kax.scatter(range(len(steps)), [0] * len(steps), s=[area(c) for c in steps], color=MUTED_INK,
                edgecolors="white", lw=0.4)
    for j, c in enumerate(steps):
        kax.text(j, -0.9, str(c), ha="center", va="top", fontsize=8, color=MUTED_INK)
    kax.text((len(steps) - 1) / 2, -2.1, f"mice with the taxon, of {k} (no dot: none)", ha="center", va="top",
             fontsize=8.5, color=MUTED_INK)
    kax.set_xlim(-0.6, len(steps) - 0.4)
    kax.set_ylim(-2.5, 0.8)
    kax.axis("off")
    scope = "taxa named no further than genus" if genus_only else "taxa of the DA section"
    fig.suptitle(f"Presence and abundance by group: the {n} {scope}", x=0.01, y=1 - 0.2 / h, ha="left",
                 fontsize=11.5, color=TEXT_INK)

    name = f"{STEM}{'_gplus' if genus_only else ''}_presence_dots"
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    absent = pres == 0
    side = {"taxa": n, "structural_zeros": int(absent.any(axis=1).sum()),
            **{f"absent_{g}": int(absent[g].sum()) for g in GROUPS},
            "min_pct": norm.vmin, "max_pct": norm.vmax}
    pd.DataFrame([side]).to_csv(os.path.join(out_dir, f"{name}.tsv"), sep="\t", index=False)
    print(f"  wrote {name}: {side}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "figures"))
    out_dir = ap.parse_args().out

    meta, pct = abundance()
    spec = pd.read_csv(META, sep="\t", index_col=0, keep_default_na=False)[["Group", "PrimaryColor"]].drop_duplicates()
    colors = dict(zip(spec.Group, spec.PrimaryColor))
    tested = pd.read_csv(os.path.join(TAB_DIR, f"{STEM}.tsv"), sep="\t").taxon
    untestable = pd.read_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t").taxon
    if set(tested) & set(untestable):
        raise SystemExit("a taxon is both tested and untestable")
    taxa = list(tested) + list(untestable)
    missing = set(taxa) - set(pct.index)
    if missing:
        raise SystemExit(f"{len(missing)} taxa not found after collapsing: {sorted(missing)[:3]}")
    n_mice = {g: int((meta.Group == g).sum()) for g in GROUPS}
    if len(set(n_mice.values())) != 1:
        raise SystemExit(f"unequal groups {n_mice}: the size key assumes one group size")
    pres = pd.DataFrame({g: (pct.loc[taxa, meta.Group == g] > 0).sum(axis=1) for g in GROUPS})
    # Mean over every mouse of the group, a mouse without the taxon counting as 0.
    mean = pd.DataFrame({g: pct.loc[taxa, meta.Group == g].mean(axis=1) for g in GROUPS})
    if (pres.loc[list(tested)] == 0).any().any() or not (pres.loc[list(untestable)] == 0).any(axis=1).all():
        raise SystemExit("tested taxa must be in every group and untestable ones missing from one")
    # One shade scale for both views, set by every taxon, so a shade means the same in each.
    shown = mean.values[pres.values > 0]
    norm = LogNorm(vmin=shown.min(), vmax=shown.max())
    for genus_only in (False, True):
        draw(genus_only, pres, mean, norm, colors, n_mice["Control"], out_dir)


if __name__ == "__main__":
    main()
