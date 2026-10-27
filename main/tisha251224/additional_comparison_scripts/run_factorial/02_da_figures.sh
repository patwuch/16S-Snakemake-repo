# Whole-design views (clustered heatmap, faceted forest) and the interaction forest,
# each also as a genus-and-above view.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "factorial views" in_sif "$QIIME" python3 "$SCRIPTS/plot_factorial_views.py"
light "factorial views, genus names" in_sif "$QIIME" python3 "$SCRIPTS/plot_factorial_views.py" --genus-names-only
light "interaction forest" in_sif "$QIIME" python3 "$SCRIPTS/plot_interaction_forest.py"
light "interaction forest, genus names" in_sif "$QIIME" python3 "$SCRIPTS/plot_interaction_forest.py" --genus-names-only
