#!/usr/bin/env bash
# ahmad260430: the Treatment x Radiation factorial report (Day 17), from the pipeline's outputs.
#
# Run the standard pipeline first. Then, from anywhere:
#   bash run_all.sh [--dry-run] [--force] [--from NN | --only NN]
#   bash run_all.sh --record-metadata
#
# tables/metadata_day17.tsv is written by hand (metadata.tsv's Day 17 rows). Every
# step stops if metadata.tsv has changed since it was recorded. Most scripts still
# read metadata.tsv and select Day 17 themselves.
# Heavy steps (count export, ALDEx3, PICRUSt2, MaAsLin3) skip while their outputs
# are newer than their inputs; every other step always reruns. PICRUSt2 is shared
# with run_day17_day25/.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh" --no-check
run_steps "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$@"
