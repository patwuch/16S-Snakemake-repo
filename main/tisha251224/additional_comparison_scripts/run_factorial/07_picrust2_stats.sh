# PICRUSt2 per-feature tests (MaAsLin3, one level and one core at a time: memory grows
# with cores), PERMANOVA, ggpicrust2 figures and NSTI.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
for level in MetaCyc KO EC; do
    OUTS=("$T/PICRUSt2_${level}_maaslin3_rows.tsv" "$T/PICRUSt2_${level}_maaslin3_rows_summary.tsv")
    INS=("$MA/metadata.tsv" "$P/pipeline" "$SCRIPTS/run_picrust2_maaslin3.R")
    heavy "MaAsLin3 $level" in_sif maaslin3-env Rscript "$SCRIPTS/run_picrust2_maaslin3.R" "$level" --cores 1
done
light "PERMANOVA" in_sif ancombc2-env Rscript "$SCRIPTS/picrust2_permanova.R"
light "ggpicrust2 PCA" in_sif ggpicrust2-env Rscript "$SCRIPTS/plot_ggpicrust2.R" pca
light "ggpicrust2 heatmaps and volcanoes" in_sif ggpicrust2-env Rscript "$SCRIPTS/plot_ggpicrust2.R" results
light "NSTI" in_sif "$QIIME" python3 "$SCRIPTS/plot_picrust2_nsti.py"
