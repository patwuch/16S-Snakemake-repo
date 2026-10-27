# Diversity statistics over both arms, then the beta-diversity figures.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "diversity tests" in_sif "$QIIME" python3 "$SCRIPTS/run_diversity_day17_day25.py" "$C"
light "beta figures" in_sif "$QIIME" python3 "$SCRIPTS/plot_beta_day17_day25.py"
