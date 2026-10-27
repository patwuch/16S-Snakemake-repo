# Structural zeros. The taxon lists (seen in every cell / absent from a cell) come
# from run_ancombc2_factorial.R's export step, which fits nothing; the per-cell
# counts from presence_interaction.R; each mouse's share in them from
# plot_untestable_heatmap.py. Then the presence dot map.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
light "presence taxon lists" in_sif ancombc2-env Rscript "$SCRIPTS/run_ancombc2_factorial.R" --export-only
light "structural zeros" in_sif "$QIIME" python3 "$SCRIPTS/structural_zeros.py"
light "presence by cell" bash -c "cd '$SCRIPTS' && apptainer exec --home '$REPO' '$REPO/containers/ancombc2-env.sif' \
    Rscript presence_interaction.R --comparison comparison-Diet-x-Condition-factorial \
    --meta comparison-Diet-x-Condition-factorial/metadata_factorial.tsv --factors Diet,Condition \
    --levels-a Cellulose,Inulin,RS --levels-b Control,IR --stem presence_factorial_Greengenes2 --firth"
light "untestable shares" in_sif "$QIIME" python3 "$SCRIPTS/plot_untestable_heatmap.py"
light "presence dot map" in_sif "$QIIME" python3 "$SCRIPTS/plot_presence_dots_factorial.py"
