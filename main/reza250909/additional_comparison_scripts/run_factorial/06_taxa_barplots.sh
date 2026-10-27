# Composition barplots by group, one per rank.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
for level in Kingdom Phylum Class Order Family Genus Species; do
    o="$F/taxa_factorial_Greengenes2_$level"
    light "taxa barplot $level" in_sif qiime2-2025.10-amplicon-core python3 "$SCRIPTS/plot_taxa_barplot.py" \
        --table-biom "$BIOM" --taxonomy-tsv "$GG2_TAX" --metadata-tsv "$META" --db Greengenes2 --level "$level" \
        --factor Group --group-order Control SCFAs IR SCFAs+IR --top-n 20 --title "Diet × Radiation, by group" \
        --output-samples-svg "${o}_samples.svg" --output-samples-png "${o}_samples.png" \
        --output-groups-svg "${o}_groups.svg" --output-groups-png "${o}_groups.png"
done
