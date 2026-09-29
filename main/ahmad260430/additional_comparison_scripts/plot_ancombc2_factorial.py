#!/usr/bin/env python3
"""Effect-size barplots and cladograms for the factorial ANCOM-BC2 coefficients.

Renders the .res files run_ancombc2_factorial.R writes, through the pipeline's
own shared renderer (scripts/effectsize_cladogram_render.py) so these charts are
directly comparable to the LEfSe ones rather than merely similar-looking.

Each .res is one two-class contrast (a main-effect coefficient conditional on
the other factor's reference level), so it gets its own pair of figures.

Colours come from the unsubsetted metadata's PrimaryColor/PrimaryOrder, keyed on
Group -- and because a collapsed feature's class is named by its Group cell, the
same palette as the diversity figures applies without a second mapping.

Usage, inside the qiime2 container:
    python3 plot_ancombc2_factorial.py [--pattern '*_minpres6_*.res']
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
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial")
# Plots live in figures/, data and tables at the root. A hundred-odd plots at
# the root bury the handful of files a person actually opens.
FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)

sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
from effectsize_cladogram_render import render_effectsize_cladogram  # noqa: E402

# group = NULL mode applies plain BH across features per coefficient. It does
# NOT use ancombc2's mdFDR machinery -- that lives in the pairwise path, which
# this analysis deliberately disables. Labelling it mdFDR would overstate it.
EFFECT_LABEL = "|log2FC| (BH-adjusted)"

TREATMENT_NAME = {"Control": "Control", "IFA": "Isoferulic acid", "Inosine": "Inosine"}
RADIATION_NAME = {"None": "unirradiated", "IR": "irradiated"}


def coef_title(stem, coef):
    """Name the comparison a coefficient actually makes in its own fit.

    Every coefficient is conditional on the fit's reference levels, which the
    file name records as `_ref<Treatment>` and `_radref<Radiation>` (Control and
    None when absent). A title built from the coefficient alone reads the same
    for "Inosine vs Control" and "Inosine vs IFA", so build it from both.
    """
    t_ref = (re.search(r"_ref(IFA|Inosine|Control)(?:_|$)", stem) or [None, "Control"])[1]
    r_ref = (re.search(r"_radref(IR|None)(?:_|$)", stem) or [None, "None"])[1]
    subset = re.search(r"_treatment(IFA|Inosine|Control)(?:_|$)", stem)
    if subset:
        t_ref = subset[1]
    # Returns (title, context). The comparison is the title; the stratum it holds
    # in goes to the kicker line above it, which keeps the title as short as the
    # old coefficient-only titles so it clears the cladogram's taxa legend.
    if coef.startswith("Treatment"):
        level = coef[len("Treatment"):]
        return f"{TREATMENT_NAME[level]} vs {TREATMENT_NAME[t_ref]}", f"{RADIATION_NAME[r_ref]} mice"
    if coef.startswith("Radiation"):
        level = coef[len("Radiation"):]
        return (f"{RADIATION_NAME[level].capitalize()} vs {RADIATION_NAME[r_ref]}",
                f"{TREATMENT_NAME[t_ref]} mice")
    return coef, None


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
    ap.add_argument("--pattern", default="ANCOMBC2_factorial_*.res")
    args = ap.parse_args()

    colors = group_palette(os.path.join(STUDY_DIR, "metadata.tsv"))
    files = sorted(glob.glob(os.path.join(TAB_DIR, args.pattern)))
    if not files:
        raise SystemExit(f"no .res matching {args.pattern} in {TAB_DIR}")

    for res_path in files:
        stem = os.path.basename(res_path)[: -len(".res")]
        m = re.match(r"ANCOMBC2_factorial_(\w+?)_(.+)_([A-Za-z]+(?:IFA|Inosine|Control|IR|None))$", stem)
        db = m.group(1) if m else "unknown"
        coef = m.group(3) if m else stem.rsplit("_", 1)[-1]
        marked = sum(1 for ln in open(res_path)
                     if len(ln.rstrip("\n").split("\t")) >= 3 and ln.split("\t")[2].strip())
        subtitle, context = coef_title(stem, coef)
        if re.search(r"_lvl6(?:_|$)", stem):
            context = f"{context} · genus level" if context else "genus level"
        if "_gplus_" in stem:
            # Same fit as the species figure; only hits without a species call are drawn.
            context = f"{context} · genus-level names only" if context else "genus-level names only"

        sidecar = res_path[: -len(".res")] + ".struczero.tsv"
        render_effectsize_cladogram(
            res_path=res_path,
            barplot_svg=os.path.join(FIG_DIR, f"{stem}_barplot.svg"),
            cladogram_svg=os.path.join(FIG_DIR, f"{stem}_cladogram.svg"),
            barplot_png=os.path.join(FIG_DIR, f"{stem}_barplot.png"),
            cladogram_png=os.path.join(FIG_DIR, f"{stem}_cladogram.png"),
            colors_map=colors,
            db=db,
            barplot_title=f"ANCOM-BC2: {subtitle}",
            cladogram_title=f"ANCOM-BC2: {subtitle}",
            effect_label=EFFECT_LABEL,
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
