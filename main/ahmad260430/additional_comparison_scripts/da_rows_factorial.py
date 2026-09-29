#!/usr/bin/env python3
"""The differential-abundance table comparison-Treatment-x-Radiation-factorial/report.html shows.

One row per question x contrast, computed from the .res files: hit count,
largest log2FC in each direction, smallest q, the taxa, and the count under the
stricter prevalence filter (--min-presence 12). The interaction rows come from
plot_interaction_forest.py's joint q instead, since interactions get no .res.
check_report_factorial.py compares the report against these values; run this
directly to print them.

Ported from tisha251224/additional_comparison_scripts/da_rows.py.

Usage:
    python3 da_rows_factorial.py
"""
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial", "tables")
PRE = "ANCOMBC2_factorial_Greengenes2_Day17"
PREV, STRICT = "_minpres6", "_minpres12"

# (question, label, first-named group, .res name after the prevalence tag)
ROWS = [
    ("Do the treatments differ in unirradiated mice?", "IFA vs Control", "IFA", "_pairwise_TreatmentIFA"),
    ("Do the treatments differ in unirradiated mice?", "Inosine vs Control", "Inosine", "_pairwise_TreatmentInosine"),
    ("Do the treatments differ in unirradiated mice?", "IFA vs Inosine", "IFA", "_pairwise_refIFA_TreatmentInosine"),
    ("Do the treatments differ in irradiated mice?", "IFA+IR vs IR", "IFA+IR", "_radrefIR_pairwise_TreatmentIFA"),
    ("Do the treatments differ in irradiated mice?", "Inosine+IR vs IR", "Inosine+IR", "_radrefIR_pairwise_TreatmentInosine"),
    ("Do the treatments differ in irradiated mice?", "IFA+IR vs Inosine+IR", "IFA+IR", "_radrefIR_pairwise_refIFA_TreatmentInosine"),
    ("Does radiation change the microbiome within each treatment?", "Control", "IR", "_RadiationIR"),
    ("Does radiation change the microbiome within each treatment?", "IFA", "IFA+IR", "_refIFA_RadiationIR"),
    ("Does radiation change the microbiome within each treatment?", "Inosine", "Inosine+IR", "_refInosine_RadiationIR"),
]


def short(lineage):
    parts = [p for p in lineage.split(".") if p and not p.endswith("__")]
    return parts[-1] if parts else lineage


def markers(path):
    out = []
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 5 and c[2].strip() and c[3].strip():
            out.append((c[0], c[2], float(c[3]), float(c[4])))
    return out


def res_path(prev, name):
    return os.path.join(TAB_DIR, f"{PRE}{prev}{name}.res")


def rows():
    for q, label, first, name in ROWS:
        m = markers(res_path(PREV, name))
        strict = markers(res_path(STRICT, name))
        up = [x for x in m if x[1] == first]
        down = [x for x in m if x[1] != first]
        yield dict(question=q, label=label, first=first, n=len(m), n_up=len(up), n_down=len(down),
                   max_up=max((x[2] for x in up), default=None),
                   max_down=max((x[2] for x in down), default=None),
                   min_q=min((x[3] for x in m), default=None),
                   n_strict=len(strict),
                   shared=len({x[0] for x in m} & {x[0] for x in strict}),
                   strict_only=sorted(short(x[0]) for x in strict if x[0] not in {y[0] for y in m}),
                   up_taxa=sorted(short(x[0]) for x in up), down_taxa=sorted(short(x[0]) for x in down))


if __name__ == "__main__":
    for r in rows():
        f = lambda v: "—" if v is None else f"{v:.2f}"
        print(f"{r['label']:<22} n={r['n']:<2} up {r['n_up']} (max +{f(r['max_up'])}) down {r['n_down']} (max −{f(r['max_down'])}) "
              f"min q={'—' if r['min_q'] is None else format(r['min_q'], '.1e')}  >=12: {r['n_strict']} (shared {r['shared']})")
        print(f"    up:   {', '.join(r['up_taxa'])}\n    down: {', '.join(r['down_taxa'])}")
        if r["strict_only"]:
            print(f"    only under the stricter filter: {', '.join(r['strict_only'])}")
