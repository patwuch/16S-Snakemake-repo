# Composition barplots by group, one per rank, drawn with reza250909's plot_taxa_barplot.py
# (ahmad's own copy appends " arm" to the title).
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
for level in Kingdom Phylum Class Order Family Genus Species; do
    o="$F/taxa_factorial_Greengenes2_$level"
    light "taxa barplot $level" in_sif "$QIIME" python3 "$REPO/main/reza250909/additional_comparison_scripts/plot_taxa_barplot.py" \
        --table-biom "$BIOM" --taxonomy-tsv "$GG2_TAX" --metadata-tsv "$META" --db Greengenes2 --level "$level" \
        --factor Group --group-order Control IFA Inosine IR IFA+IR Inosine+IR --title "Treatment × radiation, Day 17" \
        --output-samples-svg "${o}_samples.svg" --output-samples-png "${o}_samples.png" \
        --output-groups-svg "${o}_groups.svg" --output-groups-png "${o}_groups.png"
done
