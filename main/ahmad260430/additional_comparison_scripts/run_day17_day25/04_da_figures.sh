# Structural zeros (also writes metadata_joint_cells.tsv), the presence dot map, and the
# MaAsLin3 heatmaps and forest plots.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "structural zeros" in_sif "$QIIME" python3 "$SCRIPTS/structural_zeros_day17_day25.py"
light "presence dot map" in_sif "$QIIME" python3 "$SCRIPTS/plot_presence_day17_day25.py"
light "MaAsLin3 figures" in_sif "$QIIME" python3 "$SCRIPTS/plot_maaslin3_day17_day25.py"
