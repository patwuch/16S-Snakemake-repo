# Paths for every step of ahmad260430's Day 17 vs Day 25 comparison (IFA and IFA+IR arms).
# Sourced by run_all.sh and by each step; stops unless metadata.tsv still matches the
# version recorded for metadata_joint_cells.tsv (pass --no-check to skip that).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)/scripts/project_run_lib.sh"

S="$REPO/main/ahmad260430"
SCRIPTS="$S/additional_comparison_scripts"
C="$S/comparison-day17-day25"
T="$C/tables" F="$C/figures" P="$S/picrust2"
CORE_METRICS="$S/qiime2_artifacts/core-metrics-results"

STUDY_META="$S/metadata.tsv"
META="$C/metadata_joint_cells.tsv"
BIOM="$S/tables/study-seqs.biom"
FNA="$S/tables/study-seqs.fna"
GG2_TAX="$S/tables/exported-taxonomy/Greengenes2_taxonomy.tsv"
QIIME=qiime2-2025.10-amplicon-core
ARMS=(IFA IFA_IR)
declare -A ARM_GROUP=([IFA]="IFA" [IFA_IR]="IFA+IR")
DAY_COLORS=(--color Day17=#E69F00 --color Day25=#56B4E9)

[ "${1:-}" = "--no-check" ] || check_metadata "$META" "$STUDY_META"
