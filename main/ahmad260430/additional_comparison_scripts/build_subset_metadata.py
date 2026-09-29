#!/usr/bin/env python3
"""Write a metadata.tsv subset for one treatment arm's Day17 vs Day25 samples.

'Group' in the project metadata is the bare treatment name (e.g. 'IFA'), and
the separate 'Timepoint' column carries Day17/Day25 — so isolating a
Day17-vs-Day25 comparison within one treatment arm means filtering to that
arm's Group value alone. Pure stdlib — runs on the host, no container needed.
"""
import argparse
import csv
from pathlib import Path

DAY_ORDER = ["Day17", "Day25"]
TIMEPOINT_COLUMN = "Timepoint"


def build_subset(metadata_path: Path, treatment: str, output_path: Path) -> Path:
    with open(metadata_path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header, data_rows = rows[0], rows[1:]
    group_idx = header.index("Group")
    timepoint_idx = header.index(TIMEPOINT_COLUMN)

    kept = [r for r in data_rows if r[group_idx] == treatment and r[timepoint_idx] in DAY_ORDER]
    if not kept:
        raise ValueError(
            f"No samples found for treatment {treatment!r} "
            f"(looked for Group == {treatment!r} with {TIMEPOINT_COLUMN} in {DAY_ORDER}) "
            f"in {metadata_path}"
        )

    present_days = {r[timepoint_idx] for r in kept}
    missing_days = set(DAY_ORDER) - present_days
    if missing_days:
        raise ValueError(
            f"Treatment {treatment!r} is missing {sorted(missing_days)} samples in {metadata_path} "
            f"— only has {sorted(present_days)}. Not a Day17-vs-Day25 comparison yet."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(kept)

    print(f"Wrote {len(kept)} samples for {treatment!r} -> {output_path}")
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--treatment", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    build_subset(args.metadata, args.treatment, args.output)
