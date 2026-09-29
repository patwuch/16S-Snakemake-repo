#!/usr/bin/env python3
"""Summarize the MaAsLin 3 fits run_maaslin3_factorial.R writes.

For each fit and model term: the hits (individual q < 0.05, no fitting error),
their direction and size, and two checks from MaAsLin 3's own guidance:

  abundance hits  -- at least 10 non-zero samples in each compared category
  prevalence hits -- at least 10 samples with the taxon present and 10 with
                     it absent in each compared category

Abundance hits are also matched against ANCOM-BC2's factorial fit for the same
coefficient, where one exists.

Usage, inside the qiime2 core container:
    python3 summarize_maaslin3.py
"""
import os

import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
MA_DIR = os.path.join(OUT_DIR, "tables", "maaslin3")
ALPHA = 0.05


def short(lineage):
    parts = [p.strip() for p in lineage.split(";") if p.strip() and not p.strip().endswith("__")]
    return parts[-1] if parts else lineage


def ancombc2_hits(coef):
    df = pd.read_csv(os.path.join(OUT_DIR, "tables", "ANCOMBC2_factorial_Greengenes2.tsv"), sep="\t")
    if f"q_{coef}" not in df:
        return None
    ok = (df[f"q_{coef}"] < ALPHA) & df[f"passed_ss_{coef}"].astype(str).eq("True")
    return set(df.taxon[ok])


def main():
    key = pd.read_csv(os.path.join(MA_DIR, "feature_key.tsv"), sep="\t").set_index("feature").taxon
    counts = pd.read_csv(os.path.join(MA_DIR, "counts_Greengenes2.tsv"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(MA_DIR, "metadata.tsv"), sep="\t", index_col=0)
    present = counts > 0

    for fit in ("interaction", "additive"):
        res = pd.read_csv(os.path.join(MA_DIR, f"fit_{fit}", "all_results.tsv"), sep="\t")
        tested = res.feature.nunique()
        errors = res.error.notna().sum()
        print(f"\n=== {fit}: {tested} features, {errors} of {len(res)} coefficients carry a fitting error")
        res = res[res.metadata != "log10_reads"]
        for (model, name), grp in res.groupby(["model", "name"], sort=False):
            hits = grp[(grp.qval_individual < ALPHA) & grp.error.isna()]
            nominal = int((grp.pval_individual < ALPHA).sum())
            print(f"  {model:<9} {name:<32} hits {len(hits):>2}   (raw p < 0.05: {nominal} of {len(grp)})")
            for r in hits.sort_values("qval_individual").itertuples():
                # categories the coefficient compares: its own level vs the reference
                var = r.metadata.split(":")[0]
                level = str(r.value).split(":")[0]
                ref = {"Diet": "Cellulose", "Condition": "Control"}.get(var)
                cats = [c for c in (level, ref) if c in set(meta.get(var, pd.Series(dtype=str)))]
                pres = present[r.feature]
                if model == "prevalence":
                    ok = all(pres[meta[var] == c].sum() >= 10 and (~pres[meta[var] == c]).sum() >= 10 for c in cats) if cats else False
                else:
                    ok = all(pres[meta[var] == c].sum() >= 10 for c in cats) if cats else False
                flag = "meets guidance" if ok else "below guidance"
                warn = f"  [{r.error}]" if isinstance(r.error, str) else ""
                print(f"      {short(key[r.feature]):<45} coef {r.coef:+6.2f}  q {r.qval_individual:.2g}  "
                      f"present in {int(pres.sum())}/19  {flag}{warn}")
            if model == "abundance" and fit == "interaction":
                a = ancombc2_hits(name)
                if a is not None:
                    m = {key[f] for f in hits.feature}
                    print(f"      vs ANCOM-BC2 {name}: {len(m & a)} shared, {len(m - a)} MaAsLin 3 only, {len(a - m)} ANCOM-BC2 only")


if __name__ == "__main__":
    main()
