# Paths for every step of tisha251224's Diet x Condition factorial analysis.
# Sourced by run_all.sh and by each step; stops unless metadata.tsv still matches
# the version metadata_factorial.tsv was made from (pass --no-check to skip that).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)/scripts/project_run_lib.sh"

S="$REPO/main/tisha251224"
SCRIPTS="$S/additional_comparison_scripts"
C="$S/comparison-Diet-x-Condition-factorial"
T="$C/tables" F="$C/figures" P="$C/picrust2" MA="$C/tables/maaslin3"

STUDY_META="$S/metadata.tsv"
META="$C/metadata_factorial.tsv"
BIOM="$S/tables/study-seqs.biom"
FNA="$S/tables/study-seqs.fna"
GG2_TAX="$S/tables/exported-taxonomy/Greengenes2_taxonomy.tsv"
QIIME=qiime2-2025.10-amplicon-core

[ "${1:-}" = "--no-check" ] || check_metadata "$META" "$STUDY_META"
