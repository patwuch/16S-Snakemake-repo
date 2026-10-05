#!/usr/bin/env bash
# PICRUSt2 2.6.3 functional prediction for ahmad260430, from the raw ASV files.
# One run on all 32 samples serves both reports (the Treatment x Radiation
# factorial and the Day 17 to Day 25 comparison).
#
# Inputs are the pipeline's exported ASV sequences and ASV count table
# (tables/study-seqs.fna, tables/study-seqs.biom: 2,325 ASVs x 32 samples,
# unrarefied, before any taxonomic collapse). PICRUSt2 places each ASV in its
# reference trees, predicts gene content from the nearest reference genomes,
# weights by the ASV's counts and 16S copy number, and infers MetaCyc pathways
# from the predicted EC numbers.
#
# 2.6.x uses the PICRUSt2-SC reference database (separate bacterial and archaeal
# references). The pipeline's rules/PICRUSt2.smk still runs 2.5.2 with the old
# IMG database; its output, where present, is not used here.
#
# Defaults kept: ASVs with NSTI > 2 are excluded, max parsimony hidden-state
# prediction. --per_sequence_contrib is not needed for the unstratified tables.
# One change: sequences are placed with SEPP (-t sepp), not the default epa-ng.
# The SC bacterial reference tree has 26,868 genomes; epa-ng holds likelihood
# buffers for the whole tree at once and was killed after 18 minutes on this
# 28 GB host when run for reza250909 (no output, empty .jplace). SEPP, the
# other placement tool
# PICRUSt2 supports, splits the tree into subsets and needs far less memory.
# (picrust2_pipeline.py does not expose epa-ng's chunk size, which in any case
# does not bound the tree-sized buffers.)
#
# Usage (from anywhere):
#   bash run_picrust2.sh [threads]
set -euo pipefail

THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "$THIS_DIR")"
REPO_ROOT="$(dirname "$(dirname "$STUDY_DIR")")"
OUT_DIR="$STUDY_DIR/picrust2"
SIF="$REPO_ROOT/containers/picrust2-2.6-env.sif"
THREADS="${1:-12}"

if [ -e "$OUT_DIR/pipeline" ]; then
    echo "$OUT_DIR/pipeline exists; move it aside to rerun" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"
run() { apptainer exec --bind "$REPO_ROOT:$REPO_ROOT" "$SIF" "$@"; }

run picrust2_pipeline.py --version
run picrust2_pipeline.py \
    -s "$STUDY_DIR/tables/study-seqs.fna" \
    -i "$STUDY_DIR/tables/study-seqs.biom" \
    -o "$OUT_DIR/pipeline" \
    -p "$THREADS" \
    -t sepp

mkdir -p "$OUT_DIR/described"
run add_descriptions.py -i "$OUT_DIR/pipeline/pathways_out/path_abun_unstrat.tsv.gz" -m METACYC \
    -o "$OUT_DIR/described/pathway_abun_unstrat_described.tsv.gz"
run add_descriptions.py -i "$OUT_DIR/pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" -m EC \
    -o "$OUT_DIR/described/EC_metagenome_unstrat_described.tsv.gz"
run add_descriptions.py -i "$OUT_DIR/pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" -m KO \
    -o "$OUT_DIR/described/KO_metagenome_unstrat_described.tsv.gz"
echo "PICRUSt2 done: $OUT_DIR"
