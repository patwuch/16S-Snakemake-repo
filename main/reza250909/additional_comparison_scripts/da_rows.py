#!/usr/bin/env python3
"""The differential-abundance table report.html shows, computed from the .res files.

One row per question x contrast: hit count, largest log2FC in each direction,
smallest q, the taxa themselves, and the count at each stricter prevalence
filter of the sweep (--min-presence 8, 12, 16). check_report.py compares the
report against these values; run this directly to print them.

Ported from tisha251224/additional_comparison_scripts/da_rows.py.

Usage:
    python3 da_rows.py
"""
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial", "tables")
PRE = "ANCOMBC2_factorial_Greengenes2"
SWEEP = (8, 12, 16)

# (question, label, first-named group, .res name after the prevalence tag)
ROWS = [
    ("Does the diet change the microbiome in unirradiated mice?", "SCFAs vs Control", "SCFAs", "_DietSCFAs"),
    ("Does the diet change the microbiome in irradiated mice?", "SCFAs+IR vs IR", "SCFAs+IR", "_radrefIR_DietSCFAs"),
    ("Does radiation change the microbiome within each diet?", "IR vs Control", "IR", "_RadiationIR"),
    ("Does radiation change the microbiome within each diet?", "SCFAs+IR vs SCFAs", "SCFAs+IR", "_refSCFAs_RadiationIR"),
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
        m = markers(res_path("", name))
        up = [x for x in m if x[1] == first]
        down = [x for x in m if x[1] != first]
        sweep = {n: markers(res_path(f"_minpres{n}", name)) for n in SWEEP}
        yield dict(question=q, label=label, n=len(m), n_up=len(up), n_down=len(down),
                   max_up=max((x[2] for x in up), default=None),
                   max_down=max((x[2] for x in down), default=None),
                   min_q=min((x[3] for x in m), default=None),
                   sweep={n: len(s) for n, s in sweep.items()},
                   sweep_kept={n: sorted(short(x[0]) for x in s if x[0] in {y[0] for y in m}) for n, s in sweep.items()},
                   sweep_new={n: sorted(short(x[0]) for x in s if x[0] not in {y[0] for y in m}) for n, s in sweep.items()},
                   up_taxa=sorted(short(x[0]) for x in up), down_taxa=sorted(short(x[0]) for x in down))


if __name__ == "__main__":
    for r in rows():
        f = lambda v: "—" if v is None else f"{v:.2f}"
        print(f"{r['label']:<20} n={r['n']:<2} up {r['n_up']} (max +{f(r['max_up'])}) down {r['n_down']} (max −{f(r['max_down'])}) "
              f"min q={'—' if r['min_q'] is None else format(r['min_q'], '.2g')}  sweep {' / '.join(str(r['sweep'][n]) for n in SWEEP)}")
        print(f"    up:   {', '.join(r['up_taxa'])}\n    down: {', '.join(r['down_taxa'])}")
        for n in SWEEP:
            if r["sweep"][n]:
                print(f"    >= {n}: kept {', '.join(r['sweep_kept'][n]) or '-'}; new {', '.join(r['sweep_new'][n]) or '-'}")
