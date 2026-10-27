# Per arm: the arm's Day 17 and Day 25 samples, alpha diversity by day, and composition
# barplots (taxa ranked over both arms so they share taxa and colours).
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
for arm in "${ARMS[@]}"; do
    light "$arm samples" python3 "$SCRIPTS/build_subset_metadata.py" --metadata "$STUDY_META" \
        --treatment "${ARM_GROUP[$arm]}" --output "$T/metadata_subset_$arm.tsv"
done
for arm in "${ARMS[@]}"; do
    sub="$T/metadata_subset_$arm.tsv"
    for metric in shannon evenness simpson chao1 faith_pd; do
        if [ "$metric" = faith_pd ]; then db=Greengenes2 vec="$CORE_METRICS/Greengenes2-faith-pd-vector.qza" name="Greengenes2_$metric"
        else db="" vec="$CORE_METRICS/$metric-vector.qza" name="$metric"; fi
        light "$arm alpha $metric" in_sif "$QIIME" python3 "$SCRIPTS/plot_diversity.py" alpha --vector "$vec" \
            --metadata "$sub" --metric "$metric" --factor Timepoint --db "$db" --group-order Day17 Day25 "${DAY_COLORS[@]}" \
            --output-svg "$F/alpha_day17day25_${arm}_$name.svg" --output-png "$F/alpha_day17day25_${arm}_$name.png"
    done
    for level in Kingdom Phylum Class Order Family Genus Species; do
        o="$F/taxa_day17day25_${arm}_Greengenes2_$level"
        light "$arm taxa barplot $level" in_sif "$QIIME" python3 "$SCRIPTS/plot_taxa_barplot.py" \
            --table-biom "$BIOM" --taxonomy-tsv "$GG2_TAX" --metadata-tsv "$sub" --db Greengenes2 --level "$level" \
            --factor Timepoint --group-order Day17 Day25 --top-n 20 --arm "${ARM_GROUP[$arm]}" \
            --rank-metadata-tsv "$T/metadata_subset_IFA.tsv" "$T/metadata_subset_IFA_IR.tsv" \
            --output-samples-svg "${o}_samples.svg" --output-samples-png "${o}_samples.png" \
            --output-groups-svg "${o}_groups.svg" --output-groups-png "${o}_groups.png"
    done
done
