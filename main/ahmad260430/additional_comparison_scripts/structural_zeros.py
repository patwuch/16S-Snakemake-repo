#!/usr/bin/env python3
"""Presence/absence tests for the report's structural-zero taxa.

A structural zero, in ANCOM-BC2's sense, is a taxon with zero counts in every
sample of a group; verify_structural_zeros.R checks that the report's
<STEM>_untestable_absent_from_a_cell.tsv is exactly ancombc2's own list. The
package declares such a taxon differentially abundant between that group and
any group where it occurs, without a test. This adds one: for every contrast in
the report's DA table and every structural-zero taxon absent from one side of
it and seen on the other, an exact test of presence counts.

  * independent groups: Fisher's exact test, two-sided;
  * the same animals on both sides (PAIR_BY set): exact McNemar test on the
    animals present on one side only.

q is Benjamini-Hochberg over every (taxon, contrast) test in the report, one
family. Presence is recomputed from the counts and checked against the
untestable table.

Writes tables/<STEM>_structural_zeros.tsv, one row per test, and prints one
summary line per contrast, which check_report.py compares with the report.

Copied from reza250909/additional_comparison_scripts/structural_zeros.py; ahmad260430
holds two analyses, so the configuration is picked with --comparison.

Usage, inside the qiime2 container:
    python3 structural_zeros.py --comparison factorial|day17-day25
"""
import argparse
import os

import biom
import numpy as np
import pandas as pd
from scipy.stats import binom, fisher_exact
from statsmodels.stats.multitest import multipletests

from plot_untestable_heatmap import collapse

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)

# ---- Study configuration: ahmad260430 holds two analyses, picked by --comparison ----
CONFIGS = {
    "factorial": dict(
        OUT_DIR=os.path.join(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial"),
        META=os.path.join(STUDY_DIR, "metadata.tsv"), META_FILTER=("Timepoint", "Day17"),
        STEM="DA_factorial_Greengenes2_Day17_minpres6", PAIR_BY=None,
        CONTRASTS=[
            ("Do the treatments differ in unirradiated mice?", "IFA vs Control", "IFA", "Control"),
            ("Do the treatments differ in unirradiated mice?", "Inosine vs Control", "Inosine", "Control"),
            ("Do the treatments differ in unirradiated mice?", "IFA vs Inosine", "IFA", "Inosine"),
            ("Do the treatments differ in irradiated mice?", "IFA+IR vs IR", "IFA+IR", "IR"),
            ("Do the treatments differ in irradiated mice?", "Inosine+IR vs IR", "Inosine+IR", "IR"),
            ("Do the treatments differ in irradiated mice?", "IFA+IR vs Inosine+IR", "IFA+IR", "Inosine+IR"),
            ("Does radiation change the microbiome within each treatment?", "Control", "IR", "Control"),
            ("Does radiation change the microbiome within each treatment?", "IFA", "IFA+IR", "IFA"),
            ("Does radiation change the microbiome within each treatment?", "Inosine", "Inosine+IR", "Inosine"),
        ]),
    "day17-day25": dict(
        OUT_DIR=os.path.join(STUDY_DIR, "comparison-day17-day25"),
        META=os.path.join(STUDY_DIR, "comparison-day17-day25", "metadata_joint_cells.tsv"), META_FILTER=None,
        STEM="ANCOMBC2_day17day25_joint_Greengenes2", PAIR_BY="MouseID",
        CONTRASTS=[
            ("Does the community change from Day 17 to Day 25 within each arm?", "IFA", "IFA Day25", "IFA Day17"),
            ("Does the community change from Day 17 to Day 25 within each arm?", "IFA+IR", "IFA+IR Day25", "IFA+IR Day17"),
        ]),
}
OUT_DIR = META = META_FILTER = STEM = PAIR_BY = CONTRASTS = TAB_DIR = None


def configure(name):
    global OUT_DIR, META, META_FILTER, STEM, PAIR_BY, CONTRASTS, TAB_DIR
    for k, v in CONFIGS[name].items():
        globals()[k] = v
    TAB_DIR = os.path.join(OUT_DIR, "tables")
# ----



def mcnemar_exact(b, c):
    """Two-sided exact McNemar p from the two discordant counts."""
    n = b + c
    return 1.0 if n == 0 else min(1.0, 2 * binom.cdf(min(b, c), n, 0.5))


def presence():
    meta = pd.read_csv(META, sep="\t", index_col=0, keep_default_na=False, na_values=[""])
    if META_FILTER:
        meta = meta[meta[META_FILTER[0]] == META_FILTER[1]]
    counts = biom.load_table(os.path.join(STUDY_DIR, "tables", "study-seqs.biom")).to_dataframe(dense=True)[meta.index]
    tax = pd.read_csv(os.path.join(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                      sep="\t", index_col=0)["Taxon"]
    return meta, collapse(counts, tax) > 0


def tests():
    meta, pres = presence()
    un = pd.read_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t")
    cols = {c.split(" (n=")[0]: c for c in un.columns if "(n=" in c}
    for g, c in cols.items():
        got = pres.loc[un.taxon, meta.index[meta.Group == g]].sum(axis=1).values
        if not np.array_equal(got, un[c].values):
            raise SystemExit(f"presence in {g} disagrees with the untestable table")
    out = []
    for q, label, a, b in CONTRASTS:
        sa, sb = meta.index[meta.Group == a], meta.index[meta.Group == b]
        for t in un.taxon:
            ka, kb = int(pres.loc[t, sa].sum()), int(pres.loc[t, sb].sum())
            if (ka == 0) == (kb == 0):
                continue   # absent from both or present in both: not a structural-zero contrast
            if PAIR_BY:
                ma = pres.loc[t, sa].groupby(meta.loc[sa, PAIR_BY].values).any()
                mb = pres.loc[t, sb].groupby(meta.loc[sb, PAIR_BY].values).any()
                ma, mb = ma.align(mb, join="inner")
                p = mcnemar_exact(int((ma & ~mb).sum()), int((~ma & mb).sum()))
                test = "McNemar exact"
            else:
                p = fisher_exact([[ka, len(sa) - ka], [kb, len(sb) - kb]])[1]
                test = "Fisher exact"
            out.append(dict(question=q, contrast=label, taxon=t, absent_from=a if ka == 0 else b,
                            present_in=b if ka == 0 else a, n_present=max(ka, kb),
                            n_group=len(sb) if ka == 0 else len(sa), test=test, p=p))
    df = pd.DataFrame(out)
    df["q"] = multipletests(df.p, method="fdr_bh")[1] if len(df) else []
    return df


def summary(df):
    rows = []
    for q, label, a, b in CONTRASTS:
        d = df[df.contrast == label]
        rows.append(dict(question=q, label=label, first=a, second=b, n=len(d),
                         absent_first=int((d.absent_from == a).sum()), absent_second=int((d.absent_from == b).sum()),
                         min_p=d.p.min() if len(d) else None, min_q=d.q.min() if len(d) else None,
                         n_sig=int((d.q < 0.05).sum()), n_half=int((d.n_present * 2 >= d.n_group).sum())))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison", choices=sorted(CONFIGS), required=True)
    configure(ap.parse_args().comparison)
    df = tests()
    df.to_csv(os.path.join(TAB_DIR, f"{STEM}_structural_zeros.tsv"), sep="\t", index=False)
    print(f"{len(df)} tests; q < 0.05: {(df.q < 0.05).sum()}; smallest p {df.p.min():.3g}, smallest q {df.q.min():.3g}")
    for r in summary(df):
        print(f"  {r['label']:<24} {r['n']:>2} taxa (absent from {r['first']} {r['absent_first']}, from {r['second']} "
              f"{r['absent_second']}); present in >= half the other group: {r['n_half']}; "
              f"min p {r['min_p'] if r['min_p'] is None else format(r['min_p'], '.2g')}, "
              f"min q {r['min_q'] if r['min_q'] is None else format(r['min_q'], '.2g')}, q<0.05: {r['n_sig']}")


if __name__ == "__main__":
    main()
