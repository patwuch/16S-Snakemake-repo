# MaAsLin3 on the species table: the median check at both filters, then every
# report row as a direct coefficient, robust calls from leave-one-mouse-out refits.
# Then the abundance heatmaps, forest plots and presence dot maps.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

OUTS=("$MA/counts_Greengenes2.tsv" "$MA/feature_key.tsv" "$MA/metadata.tsv")
INS=("$BIOM" "$GG2_TAX" "$META" "$SCRIPTS/prep_maaslin3_input.R")
heavy "MaAsLin3 input tables" in_sif "$QIIME" Rscript "$SCRIPTS/prep_maaslin3_input.R"

for m in 3 6; do
    OUTS=("$MA/MaAsLin3_factorial_${TAG}_min${m}_median_check.tsv")
    INS=("$MA/counts_Greengenes2.tsv" "$MA/metadata.tsv" "$SCRIPTS/maaslin3_common.R" "$SCRIPTS/check_maaslin3_median.R")
    heavy "median check, >= $m mice" in_sif maaslin3-env Rscript "$SCRIPTS/check_maaslin3_median.R" --min-samples "$m"
done

# Abundance is tested against the median feature only when every row passes the check.
VERDICT="$MA/MaAsLin3_factorial_${TAG}_median_check_verdict.tsv"
[ -s "$VERDICT" ] || [ "$DRY_RUN" = 1 ] || die "missing $VERDICT"
if [ -s "$VERDICT" ]; then
    if awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "pass") c = i } NR > 1 && $c != "TRUE" { bad = 1 } END { exit !bad }' "$VERDICT"; then
        MEDIAN=() SUFFIX=""
    else
        MEDIAN=(--median-comparison) SUFFIX="_median"
    fi
fi
OUTS=("$MA/MaAsLin3_factorial_${TAG}_min3${SUFFIX:-}_rows.tsv" "$MA/MaAsLin3_factorial_${TAG}_min3${SUFFIX:-}_rows_summary.tsv")
INS=("$MA/counts_Greengenes2.tsv" "$MA/metadata.tsv" "$VERDICT" "$SCRIPTS/maaslin3_common.R" "$SCRIPTS/run_maaslin3_rows.R")
heavy "MaAsLin3 rows${SUFFIX:-}" in_sif maaslin3-env Rscript "$SCRIPTS/run_maaslin3_rows.R" --min-samples 3 "${MEDIAN[@]}"

# limma on the same abundance model, variance moderated; then both methods on shuffled labels.
OUTS=("$MA/limma_factorial_${TAG}_min3_rows.tsv" "$MA/limma_factorial_${TAG}_min3_rows_summary.tsv")
INS=("$MA/MaAsLin3_factorial_${TAG}_min3${SUFFIX:-}_rows.tsv" "$SCRIPTS/maaslin3_common.R" "$SCRIPTS/limma_common.R" "$SCRIPTS/run_limma_rows.R")
heavy "limma rows" in_sif maaslin3-env Rscript "$SCRIPTS/run_limma_rows.R"
OUTS=("$MA/calibration_shuffle_${TAG}_min3.tsv" "$MA/calibration_shuffle_${TAG}_min3_summary.tsv")
INS=("$MA/counts_Greengenes2.tsv" "$MA/metadata.tsv" "$SCRIPTS/maaslin3_common.R" "$SCRIPTS/limma_common.R" "$SCRIPTS/calibrate_shuffle.R")
heavy "label-shuffle calibration (200)" in_sif maaslin3-env Rscript "$SCRIPTS/calibrate_shuffle.R" --shuffles 200
light "MaAsLin3 heatmaps, forests and presence dot maps" in_sif "$QIIME" python3 "$SCRIPTS/plot_maaslin3_factorial.py"
