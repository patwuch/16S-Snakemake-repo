#!/usr/bin/env bash
# reza250909: the Diet x Radiation factorial report, from the pipeline's outputs.
#
# Run the standard pipeline first. Then, from anywhere:
#   bash run_all.sh [--dry-run] [--force] [--from NN | --only NN]
#   bash run_all.sh --record-metadata
#
# metadata_factorial.tsv is written by hand (build_factorial_metadata.py can draft
# it). Every step stops if metadata.tsv has changed since it was recorded.
# Heavy steps (ANCOM-BC2, ALDEx3, PICRUSt2, MaAsLin3) skip while their outputs
# are newer than their inputs; every other step always reruns.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh" --no-check
run_steps "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" "$@"
