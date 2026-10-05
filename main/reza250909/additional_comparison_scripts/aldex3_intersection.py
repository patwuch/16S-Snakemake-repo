#!/usr/bin/env python3
"""Taxa whose Diet x Radiation interaction passes in both ANCOM-BC2 and ALDEx3.

Joins the ANCOM-BC2 default fit (run_ancombc2_factorial.R) with the ALDEx3
fits (run_aldex3_factorial.R) on the interaction coefficient
DietSCFAs:RadiationIR, over the same 104 taxa.

A taxon passes both when ANCOM-BC2 q < 0.05 and ALDEx3 q < 0.05 at the main
gamma, with the same sign. ANCOM-BC2's sensitivity-check flag (passed_ss) is
carried as a column, not used as a filter. ANCOM-BC2 reports natural-log fold
changes; they are converted to log2 here, the unit ALDEx3 reports.

Also checks that the ALDEx3 q < 0.05 set at the main gamma is the same in the
seed and nsample reruns, where those files exist.

Usage, inside the qiime2 container:
    python3 aldex3_intersection.py
"""
import math
import os
import sys

import pandas as pd
from scipy.stats import spearmanr

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(STUDY_DIR))
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts", "plotters"))
from placeholder_names import RANK, suffix  # noqa: E402

TAB_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Radiation-factorial", "tables")
DB = "Greengenes2"
COEF = "DietSCFAs:RadiationIR"
ALPHA = 0.05
MAIN_GAMMA = "0.5"
GAMMAS = [MAIN_GAMMA, "0"]
STABILITY_RUNS = ["_seed7", "_ns2000"]
LN2 = math.log(2)


def label(lineage):
    """Deepest named rank, with the ranks below it the classifier did not reach."""
    parts = [p.strip() for p in str(lineage).split(";") if p.strip()]
    named = [i for i, p in enumerate(parts) if not p.endswith("__")]
    if not named:
        return str(lineage)
    last = named[-1]
    missing = [RANK.get(p[:1], p[:1]) for p in parts[last + 1:]]
    return parts[last] + suffix(missing)


def aldex(gamma, run=""):
    path = os.path.join(TAB_DIR, f"ALDEx3_factorial_{DB}_gamma{gamma}{run}.tsv")
    return pd.read_csv(path, sep="\t").set_index("taxon") if os.path.exists(path) else None


def main():
    anc = pd.read_csv(os.path.join(TAB_DIR, f"ANCOMBC2_factorial_{DB}.tsv"), sep="\t").set_index("taxon")
    out = pd.DataFrame(index=anc.index)
    out["label"] = [label(t) for t in out.index]
    out["n_zero_mice"] = aldex(MAIN_GAMMA)["n_zero_mice"].reindex(out.index)
    out["ancombc2_log2fc"] = anc[f"lfc_{COEF}"] / LN2
    out["ancombc2_q"] = anc[f"q_{COEF}"]
    out["ancombc2_passed_ss"] = anc[f"passed_ss_{COEF}"].astype(str).str.upper() == "TRUE"

    for g in GAMMAS:
        a = aldex(g)
        if a is None:
            sys.exit(f"missing ALDEx3 fit for gamma {g}; run run_aldex3_factorial.R --gamma {g}")
        if set(a.index) != set(anc.index):
            sys.exit(f"ALDEx3 gamma {g} taxa differ from the ANCOM-BC2 default fit")
        out[f"aldex3_log2fc_gamma{g}"] = a[f"est_{COEF}"]
        out[f"aldex3_q_gamma{g}"] = a[f"q_{COEF}"]

    est_a = out[f"aldex3_log2fc_gamma{MAIN_GAMMA}"]
    out["sign_agrees"] = (out["ancombc2_log2fc"] > 0) == (est_a > 0)
    out["passes_both"] = ((out["ancombc2_q"] < ALPHA) & (out[f"aldex3_q_gamma{MAIN_GAMMA}"] < ALPHA)
                          & out["sign_agrees"])
    out = out.sort_values("ancombc2_q").reset_index()
    path = os.path.join(TAB_DIR, f"DA_interaction_intersection_{DB}.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"Wrote {path}")

    main_hits = set(out.loc[out[f"aldex3_q_gamma{MAIN_GAMMA}"] < ALPHA, "taxon"])
    for run in STABILITY_RUNS:
        a = aldex(MAIN_GAMMA, run)
        if a is None:
            print(f"  stability {run}: not run")
            continue
        hits = set(a.index[a[f"q_{COEF}"] < ALPHA])
        same = hits == main_hits
        print(f"  stability {run}: {len(hits)} ALDEx3 hits, identical to the main run: {same}")
        if not same:
            sys.exit("ALDEx3 hit set changes with the Monte Carlo seed or nsample; raise nsample")

    anc_hit = out["ancombc2_q"] < ALPHA
    rho = spearmanr(out["ancombc2_log2fc"], est_a).correlation
    print(f"  taxa: {len(out)}")
    print(f"  ANCOM-BC2 q < {ALPHA}: {int(anc_hit.sum())} "
          f"(passing its sensitivity check: {int((anc_hit & out['ancombc2_passed_ss']).sum())})")
    for g in GAMMAS:
        print(f"  ALDEx3 gamma {g} q < {ALPHA}: {int((out[f'aldex3_q_gamma{g}'] < ALPHA).sum())}")
    print(f"  passes both (gamma {MAIN_GAMMA}, same sign): {int(out['passes_both'].sum())}")
    print(f"  sign agreement: {int(out['sign_agrees'].sum())} of {len(out)}; "
          f"among ANCOM-BC2 hits {int(out.loc[anc_hit, 'sign_agrees'].sum())} of {int(anc_hit.sum())}")
    print(f"  Spearman rho of the two estimates: {rho:.2f}")
    zf = out["n_zero_mice"] == 0
    for name, sel in [("in all 32 mice", zf), ("zero in 1-8 mice", (out["n_zero_mice"] > 0) & (out["n_zero_mice"] <= 8)),
                      ("zero in > 8 mice", out["n_zero_mice"] > 8)]:
        s = out[sel]
        r = spearmanr(s["ancombc2_log2fc"], s[f"aldex3_log2fc_gamma{MAIN_GAMMA}"]).correlation
        print(f"  {name}: {len(s)} taxa, rho {r:.2f}, ANCOM-BC2 q-hits {int((s['ancombc2_q'] < ALPHA).sum())}")
    off = (out.loc[zf, "ancombc2_log2fc"] - out.loc[zf, f"aldex3_log2fc_gamma{MAIN_GAMMA}"])
    print(f"  offset ANCOM-BC2 - ALDEx3 over the taxa in all mice: median {off.median():+.2f} log2, "
          f"range {off.min():+.2f} to {off.max():+.2f}")
    print(f"  zero mice among ANCOM-BC2 hits: {sorted(out.loc[anc_hit, 'n_zero_mice'].astype(int).tolist())}")


if __name__ == "__main__":
    main()
