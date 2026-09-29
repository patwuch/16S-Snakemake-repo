#!/usr/bin/env python3
"""Write the metadata the Diet x Condition factorial scripts read.

metadata.tsv is left untouched. Its header ends in a stray tab (an empty sixth
column), its level names carry "%" and spaces that R turns into awkward model
column names, and it has no colour columns. This writes a clean copy beside the
analysis with:

  Diet       Cellulose / Inulin / RS        (reference Cellulose)
  Condition  Control / IR                   (Control = not irradiated)
  Group      display label per cell, e.g. "Inulin+IR"
  DietLabel  the original diet name, for figure text
  PrimaryColor, PrimaryOrder   one colour per cell, unirradiated cells first

Timepoint is dropped: the analysis treats every sample as a replicate of its
Diet x Condition cell.

Usage:
    python3 build_factorial_metadata.py
"""
import csv
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)
OUT_DIR = os.path.join(STUDY_DIR, "comparison-Diet-x-Condition-factorial")
OUT = os.path.join(OUT_DIR, "metadata_factorial.tsv")

DIET = {"0.2% cellulose": "Cellulose", "10% inulin": "Inulin", "10% resistant starch": "RS"}
CONDITION = {"Control": "Control", "IR": "IR"}
DIET_ORDER = ["Cellulose", "Inulin", "RS"]
# Light step for unirradiated, dark step of the same hue for irradiated.
# Validated with the dataviz palette checker (all pairs, light surface): every
# check passes; the three light steps sit under 3:1 contrast on white, which
# the axis labels and legends in every figure cover.
COLOR = {("Cellulose", "Control"): "#3fbf8f", ("Inulin", "Control"): "#6aa3e8",
         ("RS", "Control"): "#f08a4b", ("Cellulose", "IR"): "#12805a",
         ("Inulin", "IR"): "#1f5bb0", ("RS", "IR"): "#b8460f"}


def main():
    with open(os.path.join(STUDY_DIR, "metadata.tsv"), newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    out = []
    for r in rows:
        diet, cond = DIET[r["Diet"].strip()], CONDITION[r["Condition"].strip()]
        order = DIET_ORDER.index(diet) + (3 if cond == "IR" else 0) + 1
        out.append({"sampleid": r["sampleid"].strip(), "Group": diet if cond == "Control" else f"{diet}+IR",
                    "Diet": diet, "Condition": cond, "DietLabel": r["Diet"].strip(),
                    "PrimaryColor": COLOR[(diet, cond)], "PrimaryOrder": order})
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {OUT}: {len(out)} samples")


if __name__ == "__main__":
    main()
