#!/usr/bin/env python3
"""Draw each DA panel again with its structural zeros added.

For every panel: take the structural zeros structural_zeros.py found for that
contrast, write a figure-only copy of the panel's .res with their nodes added
(write_structural_zero_res.R) plus a sidecar naming them, derive the "Genus and
above" copy, and render both through the shared renderer's struc_zeros option
(with label_boxes). A structural zero in every sample of the group that has it
(complete separation) is drawn at full strength, one in only some of them faded.
In the bar chart a structural zero is a hatched bar at the
panel's largest tested |log2FC| with presence counts at the tip; in the
cladogram a hollow ring with a hatched, named clade. The hatched length is a
display convention: these taxa have no log2FC.

Adapted from reza250909/additional_comparison_scripts/structural_zero_figures.py:
ahmad260430 holds two analyses that draw through different wrappers, so the
configuration is per --comparison and so is the render step.

Usage, inside the qiime2 container:
    python3 structural_zero_figures.py --comparison factorial|day17-day25
"""
import argparse
import os
import re
import subprocess
import sys

import pandas as pd

import structural_zeros as sz
from genus_plus_res import species_resolved

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)

FACT = "ANCOMBC2_factorial_Greengenes2_Day17"
D1725 = "ANCOMBC2_day17day25"
TREATMENT_COEF = re.compile(r"(.+)_((?:Treatment|Radiation)[A-Za-z]+)$")

# contrast label -> [(.res stem, fitted table[, figure context for the day17-day25 wrapper])]
PANELS = {
    "factorial": {
        "IFA vs Control": [(f"{FACT}_minpres6_pairwise_TreatmentIFA", f"{FACT}_minpres6_pairwise")],
        "Inosine vs Control": [(f"{FACT}_minpres6_pairwise_TreatmentInosine", f"{FACT}_minpres6_pairwise")],
        "IFA vs Inosine": [(f"{FACT}_minpres6_pairwise_refIFA_TreatmentInosine", f"{FACT}_minpres6_pairwise")],
        "IFA+IR vs IR": [(f"{FACT}_minpres6_radrefIR_pairwise_TreatmentIFA", f"{FACT}_minpres6_radrefIR_pairwise")],
        "Inosine+IR vs IR": [(f"{FACT}_minpres6_radrefIR_pairwise_TreatmentInosine",
                              f"{FACT}_minpres6_radrefIR_pairwise")],
        "IFA+IR vs Inosine+IR": [(f"{FACT}_minpres6_radrefIR_pairwise_refIFA_TreatmentInosine",
                                  f"{FACT}_minpres6_radrefIR_pairwise")],
        "Control": [(f"{FACT}_minpres6_RadiationIR", f"{FACT}_minpres6")],
        "IFA": [(f"{FACT}_minpres6_refIFA_RadiationIR", f"{FACT}_minpres6_refIFA")],
        "Inosine": [(f"{FACT}_minpres6_refInosine_RadiationIR", f"{FACT}_minpres6_refInosine")],
    },
    "day17-day25": {
        "IFA": [(f"{D1725}_IFA_Greengenes2", f"{D1725}_IFA_Greengenes2", "IFA arm"),
                (f"{D1725}_joint_Greengenes2_armrefIFA", f"{D1725}_joint_Greengenes2_armrefIFA",
                 "IFA arm · joint Timepoint × Arm model")],
        "IFA+IR": [(f"{D1725}_IFA_IR_Greengenes2", f"{D1725}_IFA_IR_Greengenes2", "IFA+IR arm"),
                   (f"{D1725}_joint_Greengenes2_armrefIFA_IR", f"{D1725}_joint_Greengenes2_armrefIFA_IR",
                    "IFA+IR arm · joint Timepoint × Arm model")],
    },
}


def sz_stem(stem, tag):
    m = TREATMENT_COEF.match(stem)
    return f"{m.group(1)}_{tag}_{m.group(2)}" if m else f"{stem}_{tag}"


def gplus_copy(res, side, res_g, side_g):
    """Genus-and-above view: unmark species-resolved markers and drop species-resolved rings."""
    with open(res) as fin, open(res_g, "w") as fout:
        for line in fin:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and c[2].strip() and species_resolved(c[0]):
                c = [c[0], c[1], "", "", "-"]
            fout.write("\t".join(c) + "\n")
    with open(side) as fin, open(side_g, "w") as fout:
        for line in fin:
            k = line.split("\t")[0]
            if k == "#classes" or not species_resolved(k):
                fout.write(line)


def mark_separation(side):
    """Add a 4th sidecar column: "complete" when the taxon is in every sample of the
    group that has it ("4/4 vs 0/4"), "partial" otherwise ("1/4 vs 0/4"). The
    renderer draws partial rows faded; the chance calibration in
    structural_zero_evidence.py counts the complete ones."""
    out = []
    for line in open(side):
        f = line.rstrip("\n").split("\t")
        if f[0] and f[0] != "#classes":
            k, n = f[2].split(" vs ")[0].split("/")
            f = f[:3] + ["complete" if k == n else "partial"]
        out.append("\t".join(f))
    with open(side, "w") as fh:
        fh.write("\n".join(out) + "\n")


def render_day17(stem, context):
    """Same call run_comparison.py's render_ancombc2 makes, plus the sidecar."""
    fig = os.path.join(STUDY_DIR, "comparison-day17-day25", "figures")
    tab = sz.TAB_DIR
    colors = ["--color", "Day17=#E69F00", "--color", "Day25=#56B4E9"]
    subprocess.run([sys.executable, os.path.join(THIS_DIR, "plot_lefse_res.py"),
                    "--res", os.path.join(tab, f"{stem}.res"), "--db", "Greengenes2", *colors,
                    "--lda-title", "ANCOM-BC2 Effect Size", "--cladogram-title", "ANCOM-BC2 Cladogram",
                    "--context", context, "--effect-label", "|log2FC| (BH-adjusted)",
                    "--struc-zeros", os.path.join(tab, f"{stem}.struczero.tsv"), "--label-boxes", "--resolve-placeholders",
                    "--legend-ncol", "2",
                    "--lda-svg", os.path.join(fig, f"{stem}_barplot.svg"),
                    "--cladogram-svg", os.path.join(fig, f"{stem}_cladogram.svg"),
                    "--lda-png", os.path.join(fig, f"{stem}_barplot.png"),
                    "--cladogram-png", os.path.join(fig, f"{stem}_cladogram.png")], check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison", choices=sorted(PANELS), required=True)
    comp = ap.parse_args().comparison
    sz.configure(comp)
    tab = os.path.relpath(sz.TAB_DIR, STUDY_DIR)
    tests = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
    meta = pd.read_csv(sz.META, sep="\t", index_col=0, keep_default_na=False)
    if sz.META_FILTER:
        meta = meta[meta[sz.META_FILTER[0]] == sz.META_FILTER[1]]
    size = meta.Group.value_counts()
    plan, made = [], []
    for _q, label, first, second in sz.CONTRASTS:
        for panel in PANELS[comp][label]:
            stem, fit = panel[0], panel[1]
            out = sz_stem(stem, "sz")
            made.append((out, panel[2] if len(panel) > 2 else None))
            if comp == "factorial":
                src = dict(meta="metadata.tsv", filter="Timepoint=Day17", select_col="Group",
                           select_values=f"{first},{second}", class_col="Group", first=first, second=second)
                cls = lambda g: g
            else:
                arm = first.rsplit(" ", 1)[0]
                src = dict(meta="metadata.tsv", filter="", select_col="Group", select_values=arm,
                           class_col="Timepoint", first="Day25", second="Day17")
                cls = lambda g: g.rsplit(" ", 1)[1]
            for r in tests[tests.contrast == label].itertuples():
                plan.append(dict(res_in=f"{tab}/{stem}.res", res_out=f"{tab}/{out}.res",
                                 sidecar_out=f"{tab}/{out}.struczero.tsv", fit_tsv=f"{tab}/{fit}.tsv",
                                 **src, taxon=r.taxon, present_in=cls(r.present_in),
                                 annotation=f"{r.n_present}/{r.n_group} vs 0/{size[r.absent_from]}"))
    plan_path = os.path.join(sz.TAB_DIR, f"{sz.STEM}_sz_plan.tsv")
    pd.DataFrame(plan).to_csv(plan_path, sep="\t", index=False)
    subprocess.run(["Rscript", os.path.join(THIS_DIR, "write_structural_zero_res.R"), plan_path], check=True)
    for out, context in made:
        res = os.path.join(sz.TAB_DIR, f"{out}.res")
        if not os.path.exists(res):
            continue   # no structural zeros in this contrast
        g = out.replace("_sz_", "_sz_gplus_") if "_sz_" in out else f"{out}_gplus"
        mark_separation(res[:-4] + ".struczero.tsv")
        gplus_copy(res, res[:-4] + ".struczero.tsv", os.path.join(sz.TAB_DIR, f"{g}.res"),
                   os.path.join(sz.TAB_DIR, f"{g}.struczero.tsv"))
        if comp == "day17-day25":
            render_day17(out, context)
            render_day17(g, f"{context} · genus-level names only")
    if comp == "factorial":
        subprocess.run([sys.executable, os.path.join(THIS_DIR, "plot_ancombc2_factorial.py"),
                        "--pattern", f"{FACT}_minpres6_*sz_*.res"], check=True)


if __name__ == "__main__":
    main()
