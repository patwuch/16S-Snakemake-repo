# MaAsLin3 on both arms together (~ Timepoint * Arm + reads + strata(MouseID)), Greengenes2,
# the database the report uses.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
OUTS=("$T/collapsed_counts_day17day25_Greengenes2.tsv" "$T/samples_day17day25.tsv")
INS=("$BIOM" "$GG2_TAX" "$META.source" "$SCRIPTS/export_counts_day17_day25.R")
heavy "count export" in_sif ancombc2-env Rscript "$SCRIPTS/export_counts_day17_day25.R" --db Greengenes2

OUTS=("$T/MaAsLin3_day17day25_Greengenes2.tsv" "$T/MaAsLin3_day17day25_Greengenes2_all_sets.tsv")
INS=("$T/collapsed_counts_day17day25_Greengenes2.tsv" "$T/samples_day17day25.tsv" "$SCRIPTS/run_maaslin3_day17_day25.R")
heavy "MaAsLin3" in_sif maaslin3-env Rscript "$SCRIPTS/run_maaslin3_day17_day25.R" --db Greengenes2
