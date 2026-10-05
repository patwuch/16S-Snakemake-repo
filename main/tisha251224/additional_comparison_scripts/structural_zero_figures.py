#!/usr/bin/env python3
"""Draw each DA panel again with its structural zeros added.

For every panel in PANELS: take the structural zeros structural_zeros.py found
for that contrast, write a figure-only copy of the panel's .res with their
nodes added (write_structural_zero_res.R) plus a sidecar naming them, derive the
"Genus and above" copy, and render both through the shared renderer's
struc_zeros option (with label_boxes). In the bar chart a structural zero is a
hatched bar at the panel's largest tested |log2FC| with presence counts at the
tip; in the cladogram a hollow ring with a hatched, named clade. The hatched
length is a display convention: these taxa have no log2FC.

Outputs sit beside the originals with "_sz" in the name, before the coefficient
(".._radrefIR_sz_DietSCFAs"), so the tested-only figures stay on disk.

Copied from reza250909/additional_comparison_scripts/structural_zero_figures.py; only the
configuration block differs.

Usage, inside the qiime2 container:
    python3 structural_zero_figures.py
"""
import os
import re
import subprocess
import sys

import pandas as pd

import structural_zeros as sz
from genus_plus_res import species_resolved

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)

# ---- Study configuration ----
COMP = "comparison-Diet-x-Condition-factorial"
META, FILTER, CLASS_COL = f"{COMP}/metadata_factorial.tsv", "", "Group"
PRE = "MaAsLin3_factorial_Greengenes2"
# contrast label (as structural_zeros.CONTRASTS) -> (.res stem suffix, tested-taxa table suffix);
# write_maaslin3_res.R writes both.
PANELS = {
    "Inulin vs Cellulose": ("_DietInulin", "_DietInulin_tested"),
    "RS vs Cellulose": ("_DietRS", "_DietRS_tested"),
    "RS vs Inulin": ("_refInulin_DietRS", "_refInulin_DietRS_tested"),
    "Inulin+IR vs Cellulose+IR": ("_condrefIR_DietInulin", "_condrefIR_DietInulin_tested"),
    "RS+IR vs Cellulose+IR": ("_condrefIR_DietRS", "_condrefIR_DietRS_tested"),
    "RS+IR vs Inulin+IR": ("_condrefIR_refInulin_DietRS", "_condrefIR_refInulin_DietRS_tested"),
    "Cellulose+IR vs Cellulose": ("_ConditionIR", "_ConditionIR_tested"),
    "Inulin+IR vs Inulin": ("_refInulin_ConditionIR", "_refInulin_ConditionIR_tested"),
    "RS+IR vs RS": ("_refRS_ConditionIR", "_refRS_ConditionIR_tested"),
}
COEF = re.compile(r"(.+)_((?:Diet|Condition)[A-Za-z]+)$")
# ----


def sz_stem(stem, tag):
    m = COEF.match(stem)
    return f"{m.group(1)}_{tag}_{m.group(2)}" if m else f"{stem}_{tag}"


def render(pattern):
    subprocess.run([sys.executable, os.path.join(THIS_DIR, "plot_ancombc2_factorial.py"), "--pattern", pattern],
                   check=True)


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
    group that has it ("8/8 vs 0/8"), "partial" otherwise ("3/8 vs 0/8"). The
    renderer draws partial rows with a sparse hatch, a dashed outline and a thin
    ring, in the group's own colour. Ported from ahmad260430."""
    out = []
    for line in open(side):
        f = line.rstrip("\n").split("\t")
        if f[0] and f[0] != "#classes":
            k, n = f[2].split(" vs ")[0].split("/")
            f = f[:3] + ["complete" if k == n else "partial"]
        out.append("\t".join(f))
    with open(side, "w") as fh:
        fh.write("\n".join(out) + "\n")


def main():
    tab = os.path.relpath(sz.TAB_DIR, STUDY_DIR)
    tests = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
    meta = pd.read_csv(os.path.join(STUDY_DIR, META), sep="\t", index_col=0, keep_default_na=False)
    if FILTER:
        k, v = FILTER.split("=")
        meta = meta[meta[k] == v]
    size = meta.Group.value_counts()
    plan, stems = [], []
    for _q, label, first, second in sz.CONTRASTS:
        stem, fit = PANELS[label]
        out = sz_stem(f"{PRE}{stem}", "sz")
        stems.append(out)
        for r in tests[tests.contrast == label].itertuples():
            plan.append(dict(res_in=f"{tab}/{PRE}{stem}.res", res_out=f"{tab}/{out}.res",
                             sidecar_out=f"{tab}/{out}.struczero.tsv", fit_tsv=f"{tab}/{PRE}{fit}.tsv",
                             meta=META, filter=FILTER, select_col="Group", select_values=f"{first},{second}",
                             class_col=CLASS_COL, first=first, second=second, taxon=r.taxon,
                             present_in=r.present_in,
                             annotation=f"{r.n_present}/{r.n_group} vs 0/{size[r.absent_from]}"))
    plan_path = os.path.join(sz.TAB_DIR, f"{PRE}_sz_plan.tsv")
    pd.DataFrame(plan).to_csv(plan_path, sep="\t", index=False)
    subprocess.run(["Rscript", os.path.join(THIS_DIR, "write_structural_zero_res.R"), plan_path], check=True)
    for out in stems:
        res = os.path.join(sz.TAB_DIR, f"{out}.res")
        if not os.path.exists(res):
            continue   # no structural zeros in this contrast
        g = out.replace("_sz_", "_sz_gplus_") if "_sz_" in out else f"{out}_gplus"
        mark_separation(res[:-4] + ".struczero.tsv")
        gplus_copy(res, res[:-4] + ".struczero.tsv", os.path.join(sz.TAB_DIR, f"{g}.res"),
                   os.path.join(sz.TAB_DIR, f"{g}.struczero.tsv"))
    render(f"{PRE}_*sz_*.res")


if __name__ == "__main__":
    main()
