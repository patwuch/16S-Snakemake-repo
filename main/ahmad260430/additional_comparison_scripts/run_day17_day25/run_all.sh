#!/usr/bin/env bash
# ahmad260430: the Day 17 vs Day 25 report (IFA and IFA+IR arms), from the pipeline's outputs.
# Replaces run_comparison.py for this report.
#
# Run the standard pipeline first. Then, from anywhere:
#   bash run_all.sh [--dry-run] [--force] [--from NN | --only NN]
#   bash run_all.sh --record-metadata
#
# The analysis metadata (arm subsets, MouseID pairing, metadata_joint_cells.tsv) is
# still derived by the scripts from metadata.tsv; every step stops if metadata.tsv has
# changed since it was recorded.
# Heavy steps (count export, MaAsLin3, PICRUSt2) skip while their outputs are newer than
# their inputs; every other step always reruns. PICRUSt2 is shared with run_factorial/.
# tables/MaAsLin3_day17day25_noreads_* come from an earlier version of
# run_maaslin3_day17_day25.R (before the read-depth term) and are read as they are.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh" --no-check
run_steps "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$@"
