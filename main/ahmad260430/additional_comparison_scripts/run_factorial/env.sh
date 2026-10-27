# Paths for every step of ahmad260430's Treatment x Radiation factorial (Day 17).
# Sourced by run_all.sh and by each step; stops unless metadata.tsv still matches
# the version tables/metadata_day17.tsv was made from (pass --no-check to skip that).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)/scripts/project_run_lib.sh"

S="$REPO/main/ahmad260430"
SCRIPTS="$S/additional_comparison_scripts"
C="$S/comparison-Treatment-x-Radiation-factorial"
T="$C/tables" F="$C/figures" P="$S/picrust2" MA="$C/tables/maaslin3"

STUDY_META="$S/metadata.tsv"
META="$T/metadata_day17.tsv"
BIOM="$S/tables/study-seqs.biom"
FNA="$S/tables/study-seqs.fna"
GG2_TAX="$S/tables/exported-taxonomy/Greengenes2_taxonomy.tsv"
TAG="Greengenes2_Day17"
QIIME=qiime2-2025.10-amplicon-core

[ "${1:-}" = "--no-check" ] || check_metadata "$META" "$STUDY_META"
