# ANCOM-BC2: 3 reference codings x 4 prevalence cut-offs.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

fit() {  # fit SUFFIX "FLAGS" COEF...
    local suffix="$1" flags="$2"; shift 2
    OUTS=("$T/$STEM$suffix.tsv")
    for c in "$@"; do OUTS+=("$T/$STEM${suffix}_$c.res"); done
    [ -z "$suffix" ] && OUTS+=("$T/${STEM}_untestable_absent_from_a_cell.tsv")
    INS=("$BIOM" "$GG2_TAX" "$META" "$SCRIPTS/run_ancombc2_factorial.R")
    # shellcheck disable=SC2086
    heavy "ANCOM-BC2 ${suffix:-default}" in_sif ancombc2-env Rscript "$SCRIPTS/run_ancombc2_factorial.R" $flags
}

for n in 4 8 12 16; do
    p=""; [ "$n" = 4 ] || p="_minpres$n"
    fit "$p"           "--min-presence $n"                           DietSCFAs RadiationIR
    fit "${p}_radrefIR" "--min-presence $n --radiation-reference IR" DietSCFAs RadiationNone
    fit "${p}_refSCFAs" "--min-presence $n --diet-reference SCFAs"   DietControl RadiationIR
done
