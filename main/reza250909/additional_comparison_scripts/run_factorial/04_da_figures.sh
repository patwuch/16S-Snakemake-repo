# Differential-abundance figures: heat trees, interaction forest, presence dot map.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "heat trees" in_sif metacoder-env Rscript "$SCRIPTS/plot_heat_trees_factorial.R"
light "heat trees, genus names" in_sif metacoder-env Rscript "$SCRIPTS/plot_heat_trees_factorial.R" --genus-names-only
light "interaction forest" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_interaction_forest.py"
light "interaction forest, genus names" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_interaction_forest.py" --genus-names-only
light "presence dot map" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_presence_dots_factorial.py"
