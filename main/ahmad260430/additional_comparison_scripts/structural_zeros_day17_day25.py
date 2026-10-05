#!/usr/bin/env python3
"""Structural zeros for Day17-vs-Day25: which taxa are absent from a whole group, described on raw
relative abundance. No test is run (see ANALYSIS_LOG, 2026-10-01).

Taxa: present in >= 2 samples of either arm (the 106 of MaAsLin3's set U). Those seen in every
Arm x Timepoint cell (49) are the abundance taxa; the rest (57) are absent from at least one cell,
and each is a structural zero in at least one of the four two-group contrasts.

Writes, under comparison-day17-day25/:
  tables/DA_day17day25_Greengenes2_untestable_absent_from_a_cell.tsv  per-cell presence of the 57,
      the input plot_untestable_heatmap.py draws
  tables/DA_day17day25_Greengenes2_structural_zeros.tsv  one row per contrast x taxon absent from
      every sample of one side: mice with it on the other side, complete or partial, and its mean
      share of reads where present
  metadata_joint_cells.tsv  cell label as Group, the arm's colour, cell order, MouseID

Stops unless its 49 equal the S1 taxa of run_maaslin3_day17_day25.R.

Usage, inside the qiime2 container:
    python3 structural_zeros_day17_day25.py
"""
import os
import re

import biom
import pandas as pd

from plot_untestable_heatmap import collapse

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
OUT_DIR = os.path.join(STUDY_DIR, "comparison-day17-day25")
TAB_DIR = os.path.join(OUT_DIR, "tables")
STEM = "DA_day17day25_Greengenes2"

MIN_SAMPLE_PRESENCE = 2
ARMS = {"IFA": "IFA", "IFA_IR": "IFA+IR"}   # Arm level -> Group value
DAYS = ["Day17", "Day25"]
# The report's two-group contrasts, (label, first cell, second cell).
CONTRASTS = [("Day 17 to 25 in IFA", "IFA Day17", "IFA Day25"),
             ("Day 17 to 25 in IFA+IR", "IFA+IR Day17", "IFA+IR Day25"),
             ("IFA+IR vs IFA at Day 17", "IFA Day17", "IFA+IR Day17"),
             ("IFA+IR vs IFA at Day 25", "IFA Day25", "IFA+IR Day25")]


def mouse_id(sid):
    m = re.match(r"^E\d+_(\d+)\.(\d+)$", sid) or re.match(r"^(\d+)\.(\d+)colon$", sid)
    if not m:
        raise SystemExit(f"cannot parse a mouse id out of {sid!r}")
    return f"{m.group(1)}_{m.group(2)}"


def label(taxon):
    """Deepest rank actually assigned, as the factorial report's figures name taxa."""
    r = [p.strip() for p in taxon.split(";") if p.strip() and not p.strip().endswith("__")]
    name = r[-1][3:]
    if r[-1].startswith("s__"):
        return name.replace("_", " ")
    empty_slot = taxon.rstrip().endswith("s__")
    if r[-1].startswith("g__"):
        return f"{name} ({'species unresolved' if empty_slot else 'genus only'})"
    rank = {"f": "family", "o": "order", "c": "class"}.get(r[-1][0], "higher rank")
    # A path ending in empty slots is a separate feature from one that stops at this rank.
    return f"{name} ({rank} only{', genus unresolved' if empty_slot else ''})"


def main():
    meta_all = pd.read_csv(os.path.join(STUDY_DIR, "metadata.tsv"), sep="\t", index_col=0,
                           dtype=str, keep_default_na=False)
    counts_all = biom.load_table(os.path.join(STUDY_DIR, "tables", "study-seqs.biom")).to_dataframe(dense=True)
    meta = meta_all[meta_all.Group.isin(ARMS.values()) & meta_all.index.isin(counts_all.columns)].copy()
    meta["Arm"] = meta.Group.map({v: k for k, v in ARMS.items()})
    meta["MouseID"] = [mouse_id(s) for s in meta.index]
    if not meta.groupby(["MouseID", "Timepoint"]).size().eq(1).all() or meta.MouseID.nunique() * 2 != len(meta):
        raise SystemExit("samples are not cleanly paired by mouse")

    counts = counts_all[meta.index]
    depth = counts.sum()
    counts = counts[counts.sum(axis=1) > 0]
    tax = pd.read_csv(os.path.join(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                      sep="\t", index_col=0)["Taxon"]
    col = collapse(counts, tax)
    share = col / depth * 100

    present = col > 0
    in_arm = pd.DataFrame({a: present[meta.index[meta.Arm == a]].sum(axis=1) >= MIN_SAMPLE_PRESENCE for a in ARMS})
    cell = meta.Group + " " + meta.Timepoint
    cells = [f"{g} {d}" for g in ARMS.values() for d in DAYS]
    per_cell = pd.DataFrame({c: present[meta.index[cell == c]].sum(axis=1) for c in cells})
    in_every = (per_cell > 0).all(axis=1)
    U = in_arm.any(axis=1)
    keep, absent_somewhere = U & in_every, U & ~in_every

    m3 = pd.read_csv(os.path.join(TAB_DIR, "MaAsLin3_day17day25_Greengenes2_all_sets.tsv"), sep="\t")
    if set(col.index[keep]) != set(m3.feature[m3.set == "S1"]) or set(col.index[U]) != set(m3.feature[m3.set == "U"]):
        raise SystemExit("selection differs from run_maaslin3_day17_day25.R's U / S1")

    n = {c: int((cell == c).sum()) for c in cells}
    out = pd.DataFrame({"taxon": col.index[absent_somewhere], "n_samples": present[absent_somewhere].sum(axis=1).values})
    for c in cells:
        out[f"{c} (n={n[c]})"] = per_cell.loc[absent_somewhere, c].values
    out.to_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t", index=False)

    rows = []
    for name, a, b in CONTRASTS:
        for t in col.index[U]:
            na, nb = per_cell.at[t, a], per_cell.at[t, b]
            if (na == 0) == (nb == 0):
                continue
            absent, there = (a, b) if na == 0 else (b, a)
            s = share.loc[t, meta.index[cell == there]]
            rows.append(dict(contrast=name, taxon=t, label=label(t), absent_in=absent, present_in=there,
                             mice_present=int((s > 0).sum()), of=n[there],
                             separation="complete" if (s > 0).all() else "partial",
                             mean_share_where_present=s[s > 0].mean(), mean_share_group=s.mean(),
                             max_share=s.max()))
    sz = pd.DataFrame(rows)
    if set(sz.taxon) != set(col.index[absent_somewhere]):
        raise SystemExit("the structural zeros over the four contrasts are not the taxa absent from a cell")
    sz.to_csv(os.path.join(TAB_DIR, f"{STEM}_structural_zeros.tsv"), sep="\t", index=False, float_format="%.4g")

    order = {c: i + 1 for i, c in enumerate(cells)}
    md = pd.DataFrame({"Group": cell, "PrimaryColor": meta.PrimaryColor, "PrimaryOrder": cell.map(order),
                       "MouseID": meta.MouseID}, index=meta.index)
    md.index.name = "sampleid"
    md.sort_values(["PrimaryOrder", "MouseID"]).to_csv(os.path.join(OUT_DIR, "metadata_joint_cells.tsv"), sep="\t")
    print(f"{U.sum()} taxa in U; {keep.sum()} in every cell (= MaAsLin3 S1); {absent_somewhere.sum()} absent from a cell")
    print(sz.groupby(["contrast", "separation"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
