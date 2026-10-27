# Short-chain fatty acid and formate panels: the curated acetate, formate, propionate and
# butyrate features (picrust2_scfa_lists.py) against every PICRUSt2 MaAsLin3 row.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "SCFA and formate panels" in_sif "$QIIME" python3 "$SCRIPTS/plot_picrust2_targets.py"
