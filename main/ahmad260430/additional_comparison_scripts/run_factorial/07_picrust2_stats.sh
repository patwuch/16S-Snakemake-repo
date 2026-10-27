# PICRUSt2 for this report: MaAsLin3 per level (one at a time: memory grows with cores),
# PERMANOVA, ggpicrust2 figures and NSTI, and the inosine / isoferulic acid panels.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
for level in MetaCyc KO EC; do
    OUTS=("$T/PICRUSt2_${level}_maaslin3_rows.tsv" "$T/PICRUSt2_${level}_maaslin3_rows_summary.tsv")
    INS=("$META.source" "$P/pipeline" "$SCRIPTS/run_picrust2_maaslin3.R")
    heavy "MaAsLin3 $level" in_sif maaslin3-env Rscript "$SCRIPTS/run_picrust2_maaslin3.R" factorial "$level"
done
light "PERMANOVA" in_sif ancombc2-env Rscript "$SCRIPTS/picrust2_permanova.R" factorial
light "ggpicrust2 figures and NSTI" in_sif ggpicrust2-env Rscript "$SCRIPTS/plot_picrust2.R" factorial
light "inosine and isoferulic acid panels" in_sif "$QIIME" python3 "$SCRIPTS/plot_picrust2_targets.py" factorial
