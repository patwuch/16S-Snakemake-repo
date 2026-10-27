# PICRUSt2 2.6 from the ASV sequences and counts (SEPP placement, hours).
# run_picrust2.sh refuses to overwrite an existing run; move picrust2/ aside to redo it.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
OUTS=("$P/pipeline" "$P/described")
INS=("$FNA" "$BIOM")
heavy "PICRUSt2" bash "$SCRIPTS/run_picrust2.sh" 12
