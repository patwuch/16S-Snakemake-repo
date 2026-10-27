# MaAsLin3 on the species table: every report row as a direct coefficient, robust
# calls from leave-one-mouse-out refits. Then the .res and interaction tables.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

OUTS=("$MA/counts_Greengenes2.tsv" "$MA/feature_key.tsv" "$MA/metadata.tsv")
INS=("$BIOM" "$GG2_TAX" "$META" "$SCRIPTS/prep_maaslin3_input.R")
heavy "MaAsLin3 input tables" in_sif "$QIIME" Rscript "$SCRIPTS/prep_maaslin3_input.R"

OUTS=("$MA/MaAsLin3_factorial_Greengenes2_min2_rows.tsv" "$MA/MaAsLin3_factorial_Greengenes2_min2_rows_summary.tsv")
INS=("$MA/counts_Greengenes2.tsv" "$MA/metadata.tsv" "$SCRIPTS/run_maaslin3_rows.R")
heavy "MaAsLin3 rows" in_sif maaslin3-env Rscript "$SCRIPTS/run_maaslin3_rows.R" --min-samples 2

light "res and interaction tables" in_sif maaslin3-env Rscript "$SCRIPTS/write_maaslin3_res.R"
