# Alpha and beta diversity at Day 17: factorial tests, figures, and the beta check.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "diversity tests" in_sif "$QIIME" Rscript "$SCRIPTS/run_diversity_factorial.R" --timepoint Day17 --permutations 9999
light "diversity figures" in_sif "$QIIME" python3 "$SCRIPTS/plot_diversity_factorial.py"
light "beta check" in_sif "$QIIME" python3 "$SCRIPTS/check_beta_factorial.py"
