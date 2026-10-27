# Structural zeros: presence tests and the descriptive tables the report quotes.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "structural zeros" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/structural_zeros.py"
