# Alpha and beta diversity: factorial tests, then figures.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "diversity tests" in_sif qiime2-2025.10-amplicon-core Rscript "$SCRIPTS/run_diversity_factorial.R"
light "diversity figures" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_diversity_factorial.py"
