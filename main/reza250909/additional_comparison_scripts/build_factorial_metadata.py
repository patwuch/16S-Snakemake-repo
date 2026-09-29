#!/usr/bin/env python3
"""Write the metadata the Diet x Radiation factorial scripts read.

metadata.tsv is left untouched. It encodes the 2x2 design as two YES/NO columns
(SCFAs, IR) beside a flattened Group; this writes a copy beside the analysis
with the factors named the way the models use them:

  Diet       Control / SCFAs            (reference Control)
  Radiation  None / IR                  (None = not irradiated, the reference)
  Group      the study's own cell label: Control, SCFAs, IR, SCFAs+IR
  PrimaryColor, PrimaryOrder   carried over from metadata.tsv unchanged

Every row's Group is checked against its two factors, so a mislabelled sample
stops the build instead of landing in the wrong cell.

Usage:
    python3 build_factorial_metadata.py
"""
import csv
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
OUT = os.path.join(OUT_DIR, "metadata_factorial.tsv")

DIET = {"NO": "Control", "YES": "SCFAs"}
RADIATION = {"NO": "None", "YES": "IR"}
GROUP = {("Control", "None"): "Control", ("SCFAs", "None"): "SCFAs",
         ("Control", "IR"): "IR", ("SCFAs", "IR"): "SCFAs+IR"}


def main():
    with open(os.path.join(STUDY_DIR, "metadata.tsv"), newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    out = []
    for r in rows:
        diet, rad = DIET[r["SCFAs"].strip()], RADIATION[r["IR"].strip()]
        group = r["Group"].strip()
        if group != GROUP[(diet, rad)]:
            raise SystemExit(f"{r['sampleid']}: Group {group!r} disagrees with SCFAs={r['SCFAs']}, IR={r['IR']}")
        out.append({"sampleid": r["sampleid"].strip(), "Group": group, "Diet": diet, "Radiation": rad,
                    "PrimaryColor": r["PrimaryColor"].strip(), "PrimaryOrder": r["PrimaryOrder"].strip()})
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {OUT}: {len(out)} samples")


if __name__ == "__main__":
    main()
