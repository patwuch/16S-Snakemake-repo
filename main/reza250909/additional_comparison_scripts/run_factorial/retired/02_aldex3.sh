# ALDEx3 on the interaction, with its stability reruns, then where it agrees with ANCOM-BC2.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

run() {  # run NAME "FLAGS"
    OUTS=("$T/ALDEx3_factorial_Greengenes2_$1.tsv")
    INS=("$BIOM" "$GG2_TAX" "$META" "$T/$STEM.tsv" "$SCRIPTS/run_aldex3_factorial.R")
    # shellcheck disable=SC2086
    heavy "ALDEx3 $1" in_sif ancombc2-env Rscript "$SCRIPTS/run_aldex3_factorial.R" $2
}

run gamma0.5        "--gamma 0.5"
run gamma0          "--gamma 0"
run gamma0.5_seed7  "--gamma 0.5 --seed 7"
run gamma0.5_ns2000 "--gamma 0.5 --nsample 2000"

light "intersection" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/aldex3_intersection.py"
light "concordance figures" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_aldex3_concordance.py"
light "concordance figures, genus names" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_aldex3_concordance.py" --genus-names-only
