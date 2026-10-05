#!/usr/bin/env python3
"""Effect-size barplots and cladograms for the factorial DA coefficients.

Renders the .res files write_maaslin3_res.R (and, before it,
run_ancombc2_factorial.R) writes, through the pipeline's
own shared renderer (scripts/effectsize_cladogram_render.py) so these charts are
directly comparable to the LEfSe ones rather than merely similar-looking.

Each .res is one two-class contrast (a main-effect coefficient conditional on
the other factor's reference level), so it gets its own pair of figures.

Colours come from metadata_factorial.tsv's PrimaryColor, keyed on Group, the
same palette as the diversity figures.

Usage, inside the qiime2 container:
    python3 plot_ancombc2_factorial.py [--pattern 'MaAsLin3_factorial_*.res']
"""
import argparse
import glob
import os
import re
import sys

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(STUDY_DIR))
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Condition-factorial")
# Plots live in figures/, data and tables at the root. A hundred-odd plots at
# the root bury the handful of files a person actually opens.
FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)

sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
from effectsize_cladogram_render import render_effectsize_cladogram  # noqa: E402

# Per method: the title prefix and the effect axis. MaAsLin3's abundance effect
# is a change in the taxon's share of reads, tested against zero.
METHOD = {
    "ANCOMBC2": ("ANCOM-BC2", "|log2FC| (BH-adjusted)"),
    "MaAsLin3": ("MaAsLin3", "|log2 change in share of reads| (BH-adjusted)"),
}

COND_NAME = {"Control": "unirradiated", "IR": "irradiated"}


def group(diet, cond):
    """A cell's Group label, as metadata_factorial.tsv and every legend write it."""
    return diet if cond == "Control" else f"{diet}+IR"


def coef_title(stem, coef):
    """Name the comparison a coefficient makes by its two groups.

    The title uses the same Group labels as the figure's legend and the
    report's table ("Inulin+IR vs Cellulose+IR"), so title, legend, panel label
    and table row all read alike. Every coefficient is conditional on the fit's
    reference levels, recorded in the file name as `_ref<Diet>` and
    `_condref<Condition>` (Cellulose and Control when absent). The kicker line
    carries the stratum.
    """
    d_ref = (re.search(r"_ref(Cellulose|Inulin|RS)(?:_|$)", stem) or [None, "Cellulose"])[1]
    c_ref = (re.search(r"_condref(IR|Control)(?:_|$)", stem) or [None, "Control"])[1]
    if coef.startswith("Diet"):
        level = coef[len("Diet"):]
        return f"{group(level, c_ref)} vs {group(d_ref, c_ref)}", f"{COND_NAME[c_ref]} mice"
    level = coef[len("Condition"):]
    return f"{group(d_ref, level)} vs {group(d_ref, c_ref)}", f"{d_ref} diet"


def group_palette(meta_path):
    meta = pd.read_csv(meta_path, sep="\t", index_col=0,
                       keep_default_na=False, na_values=[""])
    spec = meta[["Group", "PrimaryColor", "PrimaryOrder"]].drop_duplicates()
    if spec.Group.duplicated().any():
        raise ValueError("a Group carries more than one PrimaryColor/PrimaryOrder")
    spec = spec.sort_values("PrimaryOrder")
    return dict(zip(spec.Group, spec.PrimaryColor))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pattern", default="MaAsLin3_factorial_*.res")
    args = ap.parse_args()

    colors = group_palette(os.path.join(OUT_DIR, "metadata_factorial.tsv"))
    # _lvl files are an earlier genus-collapsed run, not part of this analysis.
    files = sorted(f for f in glob.glob(os.path.join(TAB_DIR, args.pattern)) if "_lvl" not in f)
    if not files:
        raise SystemExit(f"no .res matching {args.pattern} in {TAB_DIR}")

    for res_path in files:
        stem = os.path.basename(res_path)[: -len(".res")]
        m = re.match(r"(ANCOMBC2|MaAsLin3)_factorial_(Greengenes2)_?(.*?)_((?:Diet|Condition)[A-Za-z]+)$", stem)
        (method, effect_label), db, coef = METHOD[m.group(1)], m.group(2), m.group(4)
        marked = sum(1 for ln in open(res_path)
                     if len(ln.rstrip("\n").split("\t")) >= 3 and ln.split("\t")[2].strip())
        subtitle, context = coef_title(stem, coef)
        if "_gplus_" in stem:
            # Same fit as the species figure; only hits without a species call are drawn.
            context = f"{context} · genus-level names only"

        sidecar = res_path[: -len(".res")] + ".struczero.tsv"
        render_effectsize_cladogram(
            res_path=res_path,
            barplot_svg=os.path.join(FIG_DIR, f"{stem}_barplot.svg"),
            cladogram_svg=os.path.join(FIG_DIR, f"{stem}_cladogram.svg"),
            barplot_png=os.path.join(FIG_DIR, f"{stem}_barplot.png"),
            cladogram_png=os.path.join(FIG_DIR, f"{stem}_cladogram.png"),
            colors_map=colors,
            db=db,
            barplot_title=f"{method}: {subtitle}",
            cladogram_title=f"{method}: {subtitle}",
            effect_label=effect_label,
            context=context,
            nice_ticks=True,
            # structural_zero_figures.py writes a sidecar beside the "_sz" copies;
            # plain .res files have none and render exactly as before.
            struc_zeros=sidecar if os.path.exists(sidecar) else None,
            # Keyed letters upright in a box of their group colour, legible over
            # hatched structural-zero wedges.
            label_boxes=True,
            # Backfilled tips ("s__COE1_sp") are named by the rank actually
            # assigned ("g__COE1 (no species call)").
            resolve_placeholders=True,
            # Two legend columns keep the saved chart as wide as its bars, so it
            # stays legible when the report or its PDF scales it to the page.
            legend_ncol=2,
        )
        print(f"  {stem}: {marked} marked node(s) -> barplot + cladogram")


if __name__ == "__main__":
    main()
