#!/usr/bin/env python3
"""How much weight the report's structural zeros can carry.

A structural zero with complete separation (present in every sample of one
group, absent from every sample of the other) is the strongest presence
pattern these group sizes can produce, but its exact test cannot go below the
design's floor, so no single call can be significant. Two further questions
can still be answered.

1. Calibration, for the set: how many complete-separation calls would the
   design produce if the group labels carried no information? The labels are
   permuted (independent groups: shuffled across samples, group sizes kept;
   paired: the two timepoints swapped within each mouse, every configuration
   enumerated) and the complete-separation calls counted per contrast over
   the DA analysis's whole prevalence-filtered taxon set, tested and untestable
   taxa together, which permuting cannot change. This says how many of the
   observed calls chance alone accounts for, not which ones.

2. Depth, per call: could the absence be under-sampling? If the taxon were at
   the smallest share it reaches in any present-group sample, each absent-group
   sample of depth D would hold D x share reads on average; the Poisson
   probability of zero reads in every one of them is exp(-sum). Using the
   smallest share, not the mean, allows for the between-animal spread a single
   Poisson ignores; it is still a screen, not a test of the group difference.
   BH across the report's complete-separation calls.

Writes tables/<STEM>_structural_zeros_calibration.tsv (per contrast) and
tables/<STEM>_structural_zeros_depth.tsv (per complete-separation call).

Usage, inside the qiime2 container:
    python3 structural_zero_evidence.py --comparison factorial|day17-day25
"""
import argparse
import itertools
import os

import biom
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import structural_zeros as sz
from plot_untestable_heatmap import collapse

NPERM = 9999
SEED = 42


def load():
    meta = pd.read_csv(sz.META, sep="\t", index_col=0, keep_default_na=False, na_values=[""])
    if sz.META_FILTER:
        meta = meta[meta[sz.META_FILTER[0]] == sz.META_FILTER[1]]
    raw = biom.load_table(os.path.join(sz.STUDY_DIR, "tables", "study-seqs.biom")).to_dataframe(dense=True)[meta.index]
    tax = pd.read_csv(os.path.join(sz.STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                      sep="\t", index_col=0)["Taxon"]
    counts = collapse(raw, tax)
    depth = raw.sum()                      # every read of the sample, the denominator of a share
    return meta, counts, depth


def universe():
    """The DA analysis's prevalence-filtered taxa: the ones it tests plus the ones it cannot."""
    tested_file = f"{sz.STEM}_armrefIFA.tsv" if sz.PAIR_BY else f"{sz.STEM}.tsv"
    tested = pd.read_csv(os.path.join(sz.TAB_DIR, tested_file), sep="\t").taxon
    untest = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_untestable_absent_from_a_cell.tsv"), sep="\t").taxon
    return list(tested) + list(untest)


def complete_calls(pres, groups, contrasts):
    """Complete-separation count per contrast: taxa in every sample of one side and none of the other."""
    out = {}
    for _q, label, a, b in contrasts:
        pa, pb = pres[:, groups == a], pres[:, groups == b]
        full_a, none_a = pa.all(axis=1), ~pa.any(axis=1)
        full_b, none_b = pb.all(axis=1), ~pb.any(axis=1)
        out[label] = int(((full_a & none_b) | (full_b & none_a)).sum())
    return out


def calibration(meta, pres):
    groups = meta.Group.values
    obs = complete_calls(pres, groups, sz.CONTRASTS)
    rng = np.random.default_rng(SEED)
    if sz.PAIR_BY:
        # Swap the two timepoint labels within each mouse; enumerate every configuration.
        mice = meta[sz.PAIR_BY].unique()
        null = []
        for flips in itertools.product([False, True], repeat=len(mice)):
            g = groups.copy()
            for m, f in zip(mice, flips):
                if f:
                    idx = np.where(meta[sz.PAIR_BY].values == m)[0]
                    g[idx] = g[idx[::-1]]
            null.append(complete_calls(pres, g, sz.CONTRASTS))
        scheme = f"all {len(null)} within-mouse swaps"
    else:
        null = [complete_calls(pres, rng.permutation(groups), sz.CONTRASTS) for _ in range(NPERM)]
        scheme = f"{NPERM} shuffles of the group labels"
    null = pd.DataFrame(null)
    rows = []
    for _q, label, a, b in sz.CONTRASTS:
        n = null[label]
        rows.append(dict(contrast=label, observed=obs[label], expected=n.mean(),
                         p_perm=((n >= obs[label]).sum() + (0 if sz.PAIR_BY else 1)) / (len(n) + (0 if sz.PAIR_BY else 1)),
                         scheme=scheme))
    tot_obs, tot_null = sum(obs.values()), null.sum(axis=1)
    rows.append(dict(contrast="all contrasts", observed=tot_obs, expected=tot_null.mean(),
                     p_perm=((tot_null >= tot_obs).sum() + (0 if sz.PAIR_BY else 1)) / (len(tot_null) + (0 if sz.PAIR_BY else 1)),
                     scheme=scheme))
    return pd.DataFrame(rows)


def depth_check(meta, counts, depth):
    calls = pd.read_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros.tsv"), sep="\t")
    calls = calls[calls.n_present == calls.n_group].copy()
    rows = []
    for r in calls.itertuples():
        sp = meta.index[meta.Group == r.present_in]
        sa = meta.index[meta.Group == r.absent_from]
        share = counts.loc[r.taxon, sp] / depth[sp]
        if (counts.loc[r.taxon, sa] > 0).any() or (counts.loc[r.taxon, sp] == 0).any():
            raise SystemExit(f"{r.taxon}: not a complete-separation call in the counts")
        expected = float((depth[sa] * share.min()).sum())
        rows.append(dict(question=r.question, contrast=r.contrast, taxon=r.taxon, present_in=r.present_in,
                         absent_from=r.absent_from, n_present=r.n_present, n_group=r.n_group,
                         mean_share_pct=share.mean() * 100, min_share_pct=share.min() * 100,
                         expected_reads_absent=expected, p_depth=float(np.exp(-expected)), p_presence=r.p))
    df = pd.DataFrame(rows)
    df["q_depth"] = multipletests(df.p_depth, method="fdr_bh")[1] if len(df) else []
    df["absence"] = np.where(df.q_depth < 0.05, "not explained by depth", "consistent with under-sampling")
    return df.sort_values(["contrast", "mean_share_pct"], ascending=[True, False])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison", choices=sorted(sz.CONFIGS), required=True)
    sz.configure(ap.parse_args().comparison)
    meta, counts, depth = load()
    taxa = universe()
    missing = set(taxa) - set(counts.index)
    if missing:
        raise SystemExit(f"{len(missing)} DA taxa not in the collapsed table, e.g. {sorted(missing)[0]}")
    pres = (counts.loc[taxa] > 0).values

    cal = calibration(meta, pres)
    cal.to_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros_calibration.tsv"), sep="\t", index=False)
    print(f"{len(taxa)} taxa; calibration by {cal.scheme.iloc[0]}")
    print(cal.drop(columns="scheme").to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    dep = depth_check(meta, counts, depth)
    dep.to_csv(os.path.join(sz.TAB_DIR, f"{sz.STEM}_structural_zeros_depth.tsv"), sep="\t", index=False)
    print(f"\n{len(dep)} complete-separation calls, {dep.taxon.nunique()} taxa; "
          f"absence not explained by depth: {(dep.q_depth < 0.05).sum()} calls, "
          f"{dep.loc[dep.q_depth < 0.05, 'taxon'].nunique()} taxa")
    short = lambda t: t.split(";")[-1] if not t.rstrip().endswith("__") else t.split(";")[-2] + " (no sp)"
    for r in dep.itertuples():
        print(f"  {r.contrast:<22} {short(r.taxon):<40} in {r.present_in:<12} {r.n_present}/{r.n_group} "
              f"mean {r.mean_share_pct:6.3f}% min {r.min_share_pct:6.3f}% exp {r.expected_reads_absent:8.1f}  {r.absence}")


if __name__ == "__main__":
    main()
