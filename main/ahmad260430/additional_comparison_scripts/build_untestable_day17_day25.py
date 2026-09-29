#!/usr/bin/env python3
"""List the taxa the Day17-vs-Day25 joint ANCOM-BC2 fit could not test.

run_ancombc2_day17_day25_joint.R tests a taxon only if it passes the per-arm
presence rule (>= 2 samples) in either arm AND is seen at least once in every
Arm x Timepoint cell. This rebuilds that selection and writes the taxa that pass
the first rule but fail the second, with per-cell presence, in the same format
the factorial analyses use, so plot_untestable_heatmap.py can draw them.

Writes, under comparison-day17-day25/:
  tables/ANCOMBC2_day17day25_joint_Greengenes2_untestable_absent_from_a_cell.tsv
  metadata_joint_cells.tsv   one row per sample: cell label as Group, the arm's
                             PrimaryColor, cell order, MouseID; mice in the same
                             order in both timepoints so paired columns line up

Stops unless the taxa it keeps are exactly the ones in the joint fit's output.

Usage, inside the qiime2 container:
    python3 build_untestable_day17_day25.py
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
STEM = "ANCOMBC2_day17day25_joint_Greengenes2"

# Same constants as run_ancombc2_day17_day25_joint.R.
MIN_SAMPLE_PRESENCE = 2
ARMS = {"IFA": "IFA", "IFA_IR": "IFA+IR"}   # Arm level -> Group value
DAYS = ["Day17", "Day25"]


def mouse_id(sid):
    m = re.match(r"^E\d+_(\d+)\.(\d+)$", sid) or re.match(r"^(\d+)\.(\d+)colon$", sid)
    if not m:
        raise SystemExit(f"cannot parse a mouse id out of {sid!r}")
    return f"{m.group(1)}_{m.group(2)}"


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
    counts = counts[counts.sum(axis=1) > 0]
    tax = pd.read_csv(os.path.join(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                      sep="\t", index_col=0)["Taxon"]
    col = collapse(counts, tax)

    present = col > 0
    in_arm = pd.DataFrame({a: present[meta.index[meta.Arm == a]].sum(axis=1) >= MIN_SAMPLE_PRESENCE for a in ARMS})
    cell = meta.Group + " " + meta.Timepoint
    cells = [f"{g} {d}" for g in ARMS.values() for d in DAYS]
    per_cell = pd.DataFrame({c: present[meta.index[cell == c]].sum(axis=1) for c in cells})
    in_every = (per_cell > 0).all(axis=1)
    keep = in_arm.any(axis=1) & in_every
    untestable = in_arm.any(axis=1) & ~in_every

    fitted = set(pd.read_csv(os.path.join(TAB_DIR, f"{STEM}_armrefIFA.tsv"), sep="\t").taxon)
    if set(col.index[keep]) != fitted:
        raise SystemExit(f"rebuilt selection ({keep.sum()} taxa) differs from the joint fit ({len(fitted)})")

    n = {c: int((cell == c).sum()) for c in cells}
    out = pd.DataFrame({"taxon": col.index[untestable], "n_samples": present[untestable].sum(axis=1).values})
    for c in cells:
        out[f"{c} (n={n[c]})"] = per_cell.loc[untestable, c].values
    out.to_csv(os.path.join(TAB_DIR, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t", index=False)

    order = {c: i + 1 for i, c in enumerate(cells)}
    md = pd.DataFrame({"Group": cell, "PrimaryColor": meta.PrimaryColor, "PrimaryOrder": cell.map(order),
                       "MouseID": meta.MouseID}, index=meta.index)
    md.index.name = "sampleid"
    md.sort_values(["PrimaryOrder", "MouseID"]).to_csv(os.path.join(OUT_DIR, "metadata_joint_cells.tsv"), sep="\t")
    print(f"{in_arm.any(axis=1).sum()} taxa pass the per-arm rule; {keep.sum()} tested (matches the joint fit); "
          f"{untestable.sum()} untestable written")


if __name__ == "__main__":
    main()
