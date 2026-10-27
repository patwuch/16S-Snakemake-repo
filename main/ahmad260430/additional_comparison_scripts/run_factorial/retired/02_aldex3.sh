# ALDEx3 on the factorial: scale gamma 0 / 0.25 / 0.5 / 1 at prevalence 6 and 12, and two
# reruns of the main fit (seed 7, 2000 Monte Carlo samples). Then its figures.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

run() {  # run MINPRES GAMMA SUFFIX EXTRA...
    local mp="$1" g="$2" sfx="$3"; shift 3
    OUTS=("$T/ALDEx3_factorial_${TAG}_minpres${mp}_gamma$g$sfx.tsv" "$T/ALDEx3_factorial_${TAG}_minpres${mp}_gamma$g${sfx}_groupmeans.tsv")
    INS=("$T/tested_counts_${TAG}_minpres$mp.tsv" "$SCRIPTS/run_aldex3_factorial.R")
    heavy "ALDEx3 prevalence $mp gamma $g$sfx" in_sif ancombc2-env Rscript "$SCRIPTS/run_aldex3_factorial.R" \
        --min-presence "$mp" --gamma "$g" "$@"
}
for g in 0 0.25 0.5 1; do
    run 6 "$g" ""
    run 12 "$g" ""
done
run 6 0.5 _seed7 --seed 7
run 6 0.5 _ns2000 --nsample 2000

light "ALDEx3 forests, heatmaps and presence heatmaps" in_sif "$QIIME" python3 "$SCRIPTS/plot_aldex3_factorial.py"
