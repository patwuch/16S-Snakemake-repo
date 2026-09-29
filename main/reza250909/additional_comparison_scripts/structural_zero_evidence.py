#!/usr/bin/env python3
"""How much weight the report's structural zeros can carry, as a set.

A structural zero is a taxon absent from every mouse of one group of a contrast
and seen in the other. Each call has its own exact test (structural_zeros.py),
but a set of calls can be weighed against chance: shuffling the group labels
9,999 times (group sizes kept) and counting the calls each shuffle produces,
over the DA analysis's whole prevalence-filtered taxon set (tested and
untestable taxa together, which shuffling cannot change). This says how many of
the observed calls chance accounts for, not which ones.

Two counts per contrast:
  * all presence calls: the taxon in at least one mouse of one group and in
    none of the other;
  * complete separations: in every mouse of one group and in none of the other.

Ported from ahmad260430/additional_comparison_scripts/structural_zero_evidence.py,
without its Poisson depth check (withdrawn there: it cannot fail at these depths).

Writes tables/<STEM>_structural_zeros_calibration.tsv.

Usage, inside the qiime2 container:
    python3 structural_zero_evidence.py
"""
import os

import biom
import numpy as np
import pandas as pd

import structural_zeros as sz
from plot_untestable_heatmap import collapse

NPERM = 9999
SEED = 42


def load():
    meta = pd.read_csv(sz.META, sep="\t", index_col=0, keep_default_na=False, na_values=[""])
    raw = biom.load_table(os.path.join(sz.STUDY_DIR, "tables", "study-seqs.biom")).to_dataframe(dense=True)[meta.index]
    tax = pd.read_csv(os.path.join(sz.STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                      sep="\t", index_col=0)["Taxon"]
    return meta, collapse(raw, tax)


def universe():
    tested = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}.tsv"), sep="\t").taxon
    untest = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_untestable_absent_from_a_cell.tsv"), sep="\t").taxon
    return list(tested) + list(untest)


def calls(pres, groups):
    """Per contrast: (all presence calls, complete separations)."""
    out = {}
    for _q, label, a, b in sz.CONTRASTS:
        pa, pb = pres[:, groups == a], pres[:, groups == b]
        any_a, any_b = pa.any(axis=1), pb.any(axis=1)
        out[label] = (int((any_a ^ any_b).sum()),
                      int(((pa.all(axis=1) & ~any_b) | (pb.all(axis=1) & ~any_a)).sum()))
    return out


def main():
    meta, counts = load()
    taxa = universe()
    missing = set(taxa) - set(counts.index)
    if missing:
        raise SystemExit(f"{len(missing)} DA taxa not in the collapsed table, e.g. {sorted(missing)[0]}")
    pres = (counts.loc[taxa] > 0).values
    groups = meta.Group.values
    obs = calls(pres, groups)
    # The observed presence calls must be exactly the report's structural-zero tests.
    tests = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
    for label, (n_all, _n_full) in obs.items():
        if n_all != int((tests.contrast == label).sum()):
            raise SystemExit(f"{label}: {n_all} presence calls, but {int((tests.contrast == label).sum())} tests")
    rng = np.random.default_rng(SEED)
    null = [calls(pres, rng.permutation(groups)) for _ in range(NPERM)]
    rows = []
    for kind, k in (("all presence calls", 0), ("complete separations", 1)):
        per = {label: np.array([n[label][k] for n in null]) for label in obs}
        for label in obs:
            o = obs[label][k]
            rows.append(dict(kind=kind, contrast=label, observed=o, expected=per[label].mean(),
                             p_perm=((per[label] >= o).sum() + 1) / (NPERM + 1)))
        tot = sum(per.values()); o = sum(obs[l][k] for l in obs)
        rows.append(dict(kind=kind, contrast="all contrasts", observed=o, expected=tot.mean(),
                         p_perm=((tot >= o).sum() + 1) / (NPERM + 1)))
    cal = pd.DataFrame(rows)
    cal.to_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros_calibration.tsv"), sep="\t", index=False)
    print(f"{len(taxa)} taxa; {NPERM} shuffles of the group labels")
    print(cal.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
