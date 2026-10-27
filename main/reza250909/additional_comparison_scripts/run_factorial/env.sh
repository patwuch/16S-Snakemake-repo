# Paths for every step of reza250909's Diet x Radiation factorial analysis.
# Sourced by run_all.sh and by each step; stops unless metadata.tsv still matches
# the version metadata_factorial.tsv was made from (pass --no-check to skip that).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)/scripts/project_run_lib.sh"

S="$REPO/main/reza250909"
SCRIPTS="$S/additional_comparison_scripts"
C="$S/comparison-Diet-x-Radiation-factorial"
T="$C/tables" F="$C/figures" P="$C/picrust2" MA="$C/tables/maaslin3"

STUDY_META="$S/metadata.tsv"
META="$C/metadata_factorial.tsv"
BIOM="$S/tables/study-seqs.biom"
FNA="$S/tables/study-seqs.fna"
GG2_TAX="$S/tables/exported-taxonomy/Greengenes2_taxonomy.tsv"
STEM="ANCOMBC2_factorial_Greengenes2"

[ "${1:-}" = "--no-check" ] || check_metadata "$META" "$STUDY_META"
