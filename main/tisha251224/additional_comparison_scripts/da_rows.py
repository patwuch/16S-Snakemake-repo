#!/usr/bin/env python3
"""The differential-abundance table report.html shows, computed from the .res files.

One row per question x contrast: hit count, largest log2FC in each direction,
smallest q, the taxa themselves, and the same count under the stricter
prevalence filter (--min-presence 12). check_report.py compares the report
against these values; run this directly to print them.

Usage:
    python3 da_rows.py
"""
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial", "tables")
PRE = "ANCOMBC2_factorial_Greengenes2"

# (question, label, first-named group, .res name after the prevalence tag)
ROWS = [
    ("Do the diets differ in unirradiated mice?", "Inulin vs Cellulose", "Inulin", "_pairwise_DietInulin"),
    ("Do the diets differ in unirradiated mice?", "RS vs Cellulose", "RS", "_pairwise_DietRS"),
    ("Do the diets differ in unirradiated mice?", "RS vs Inulin", "RS", "_pairwise_refInulin_DietRS"),
    ("Do the diets differ in irradiated mice?", "Inulin+IR vs Cellulose+IR", "Inulin+IR", "_condrefIR_pairwise_DietInulin"),
    ("Do the diets differ in irradiated mice?", "RS+IR vs Cellulose+IR", "RS+IR", "_condrefIR_pairwise_DietRS"),
    ("Do the diets differ in irradiated mice?", "RS+IR vs Inulin+IR", "RS+IR", "_condrefIR_pairwise_refInulin_DietRS"),
    ("Does irradiation change the microbiome within each diet?", "Cellulose+IR vs Cellulose", "Cellulose+IR", "_ConditionIR"),
    ("Does irradiation change the microbiome within each diet?", "Inulin+IR vs Inulin", "Inulin+IR", "_refInulin_ConditionIR"),
    ("Does irradiation change the microbiome within each diet?", "RS+IR vs RS", "RS+IR", "_refRS_ConditionIR"),
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
    # A .res from a --pairwise run carries _pairwise after any reference tag
    # of the refit, so the prevalence tag always goes first.
    return os.path.join(TAB_DIR, f"{PRE}{prev}{name}.res")


def rows():
    for q, label, first, name in ROWS:
        m = markers(res_path("", name))
        strict = markers(res_path("_minpres12", name))
        up = [x for x in m if x[1] == first]
        down = [x for x in m if x[1] != first]
        yield dict(question=q, label=label, n=len(m), n_up=len(up), n_down=len(down),
                   max_up=max((x[2] for x in up), default=None),
                   max_down=max((x[2] for x in down), default=None),
                   min_q=min((x[3] for x in m), default=None),
                   n_strict=len(strict), strict_only=sorted(short(x[0]) for x in strict if x[0] not in {y[0] for y in m}),
                   shared=len({x[0] for x in m} & {x[0] for x in strict}),
                   up_taxa=sorted(short(x[0]) for x in up), down_taxa=sorted(short(x[0]) for x in down))


if __name__ == "__main__":
    for r in rows():
        f = lambda v: "—" if v is None else f"{v:.2f}"
        print(f"{r['label']:<26} n={r['n']:<2} up {r['n_up']} (max +{f(r['max_up'])}) down {r['n_down']} (max −{f(r['max_down'])}) "
              f"min q={'—' if r['min_q'] is None else format(r['min_q'], '.2g')}  strict n={r['n_strict']} shared={r['shared']}")
        print(f"    up:   {', '.join(r['up_taxa'])}\n    down: {', '.join(r['down_taxa'])}")
        if r["strict_only"]:
            print(f"    only under the stricter filter: {', '.join(r['strict_only'])}")
