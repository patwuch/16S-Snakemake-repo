# The species table at Day 17, filtered at prevalence 6 and 12, and the taxon lists
# (tested / absent from a cell). run_ancombc2_factorial.R --export-only fits nothing.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
for mp in 6 12; do
    OUTS=("$T/tested_counts_${TAG}_minpres$mp.tsv" "$T/collapsed_counts_${TAG}_minpres$mp.tsv" "$T/DA_factorial_${TAG}_minpres$mp.tsv"
          "$T/DA_factorial_${TAG}_minpres${mp}_untestable_absent_from_a_cell.tsv")
    INS=("$BIOM" "$GG2_TAX" "$META.source" "$SCRIPTS/run_ancombc2_factorial.R")
    heavy "counts and taxon lists, prevalence $mp" in_sif ancombc2-env Rscript "$SCRIPTS/run_ancombc2_factorial.R" \
        --db Greengenes2 --min-presence "$mp" --export-only
done
