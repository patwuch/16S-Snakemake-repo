#!/usr/bin/env python3
"""Write a "genus and above" copy of each factorial ANCOM-BC2 .res file.

Same fit, same tests, same q-values: the only change is that hits resolved to a
species are unmarked, so the barplot and cladogram show just the hits whose
deepest name is a genus or higher. This is a view of the species-level results,
not a re-fit at genus level, so it never tells a second story.

A hit counts as genus-or-above when its node stops above species, or when its
species slot is the placeholder write_res.R backfills for an empty `s__`:
`s__<Genus>_sp` under a named genus, or `s__<Family>_sp` under the
`g__<Family>_unresolved` it backfills for an empty genus. Exactly those forms --
`s__UBA3263_sp001689615` is a real species.

Output sits beside the input with `_gplus` inserted before the coefficient, so
plot_ancombc2_factorial.py picks it up with its usual file-name parsing. A
stem that ends in no factorial coefficient (the Day17-vs-Day25 files) gets
`_gplus` appended instead.

Usage:
    python3 genus_plus_res.py [--dir comparison-day17-day25] [--pattern '*Greengenes2_Day17_minpres6_*.res']
"""
import argparse
import glob
import os
import re

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STUDY_DIR = os.path.dirname(THIS_DIR)


def species_resolved(node):
    last = node.split(".")[-1]
    if not last.startswith("s__"):
        return False
    parts = node.split(".")
    genus = next((p[3:] for p in reversed(parts) if p.startswith("g__")), None)
    if genus and genus.endswith("_unresolved"):
        genus = genus[: -len("_unresolved")]
    return not (genus and last == f"s__{genus}_sp")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default="comparison-Treatment-x-Radiation-factorial")
    ap.add_argument("--pattern", default="*Greengenes2_Day17_minpres6_*.res")
    args = ap.parse_args()
    tab_dir = os.path.join(STUDY_DIR, args.dir, "tables")
    for path in sorted(glob.glob(os.path.join(tab_dir, args.pattern))):
        stem = os.path.basename(path)[:-len(".res")]
        if "_gplus" in stem or "_lvl" in stem:
            continue
        m = re.match(r"(.+)_((?:Treatment|Radiation)[A-Za-z]+)$", stem)
        out = os.path.join(tab_dir, f"{m.group(1)}_gplus_{m.group(2)}.res" if m else f"{stem}_gplus.res")
        kept = dropped = 0
        with open(path) as fin, open(out, "w") as fout:
            for line in fin:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 3 and cols[2].strip():
                    if species_resolved(cols[0]):
                        cols = [cols[0], cols[1], "", "", "-"]
                        dropped += 1
                    else:
                        kept += 1
                fout.write("\t".join(cols) + "\n")
        print(f"  {os.path.basename(out)}: {kept} kept, {dropped} species-level unmarked")


if __name__ == "__main__":
    main()
