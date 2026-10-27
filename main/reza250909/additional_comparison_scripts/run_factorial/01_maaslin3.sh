# MaAsLin3 on the species table: the median check at both filters, then every
# report row as a direct coefficient, robust calls from leave-one-mouse-out refits.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
QIIME=qiime2-2025.10-amplicon-core

OUTS=("$MA/counts_Greengenes2.tsv" "$MA/feature_key.tsv" "$MA/metadata.tsv")
INS=("$BIOM" "$GG2_TAX" "$META" "$SCRIPTS/prep_maaslin3_input.R")
heavy "MaAsLin3 input tables" in_sif "$QIIME" Rscript "$SCRIPTS/prep_maaslin3_input.R"

for m in 4 8; do
    OUTS=("$MA/MaAsLin3_factorial_Greengenes2_min${m}_median_check.tsv")
    INS=("$MA/counts_Greengenes2.tsv" "$MA/metadata.tsv" "$SCRIPTS/maaslin3_common.R" "$SCRIPTS/check_maaslin3_median.R")
    heavy "median check, >= $m mice" in_sif maaslin3-env Rscript "$SCRIPTS/check_maaslin3_median.R" --min-samples "$m"
done

# Abundance is tested against the median feature only when every row passes the check.
VERDICT="$MA/MaAsLin3_factorial_Greengenes2_median_check_verdict.tsv"
[ -s "$VERDICT" ] || [ "$DRY_RUN" = 1 ] || die "missing $VERDICT"
if [ -s "$VERDICT" ]; then
    if awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "pass") c = i } NR > 1 && $c != "TRUE" { bad = 1 } END { exit !bad }' "$VERDICT"; then
        MEDIAN=() SUFFIX=""
    else
        MEDIAN=(--median-comparison) SUFFIX="_median"
    fi
fi
OUTS=("$MA/MaAsLin3_factorial_Greengenes2_min4${SUFFIX:-}_rows.tsv" "$MA/MaAsLin3_factorial_Greengenes2_min4${SUFFIX:-}_rows_summary.tsv")
INS=("$MA/counts_Greengenes2.tsv" "$MA/metadata.tsv" "$VERDICT" "$SCRIPTS/maaslin3_common.R" "$SCRIPTS/run_maaslin3_rows.R")
heavy "MaAsLin3 rows${SUFFIX:-}" in_sif maaslin3-env Rscript "$SCRIPTS/run_maaslin3_rows.R" --min-samples 4 "${MEDIAN[@]}"
