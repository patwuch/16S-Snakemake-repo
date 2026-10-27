# Structural zeros: per-contrast presence tables, then their read-depth calibration.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "structural zeros" in_sif "$QIIME" python3 "$SCRIPTS/structural_zeros.py" --comparison factorial
light "structural-zero depth" in_sif "$QIIME" python3 "$SCRIPTS/structural_zero_evidence.py" --comparison factorial
