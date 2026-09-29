# silva138_AB_V3-V4_classifier.qza was trained under QIIME2 2021.4 (scikit-learn
# 0.24.1) with extract-reads primers matching NGS.smk's cutadapt primers exactly.
# The current QIIME_CONTAINER runs scikit-learn 1.4.2, and classify-sklearn hard-
# refuses to unpickle a classifier trained on a different scikit-learn release (it
# raises rather than risk silent corruption). No official SILVA 138 classifier
# matches this exact V3V4 primer set at 1.4.2, and retraining would mean losing the
# byte-identical, already-validated model, so this rule pins a legacy container
# that reproduces the original training environment instead.
SILVA138_LEGACY_CLASSIFY_CONTAINER = str(WORKFLOW_DIR / "containers/qiime2-2021.4-silva138-classify.sif")

if region == "region_V3V4":
    # rep-seqs-analysis.qza is written by filter_seqs_to_table in the modern
    # QIIME_CONTAINER (2025.10). classify-sklearn in the pinned 2021.4 container
    # rejects that archive outright ("... is not a QIIME 2 Artifact"), so re-wrap
    # the sequences into a 2021.4-native artifact first. A FeatureData[Sequence]
    # .qza is just a zip holding <uuid>/data/dna-sequences.fasta; pull that out
    # with stdlib zipfile (writer-version agnostic) and re-import it under the
    # legacy framework. Only the reads need this round trip -- the classifier was
    # trained under 2021.4 and already loads natively.
    rule Silva138_reimport_rep_seqs_legacy:
        input:
            rep_seqs = QIIME_DIR / "rep-seqs-analysis.qza"
        output:
            rep_seqs_legacy = QIIME_DIR / "rep-seqs-analysis-legacy.qza"
        params:
            fasta = str(QIIME_DIR / "rep-seqs-analysis-legacy.fasta")
        container:
            SILVA138_LEGACY_CLASSIFY_CONTAINER
        resources:
            mem_mb = 2000
        message:
            "[Silva138] Re-wrapping analysis rep-seqs into a QIIME 2 2021.4-native artifact for legacy classify-sklearn..."
        shell:
            """
            set -euo pipefail
            python -c "import zipfile, sys; z = zipfile.ZipFile(sys.argv[1]); name = next(n for n in z.namelist() if n.endswith('/data/dna-sequences.fasta')); open(sys.argv[2], 'wb').write(z.read(name))" {input.rep_seqs} {params.fasta}
            qiime tools import \
                --type 'FeatureData[Sequence]' \
                --input-path {params.fasta} \
                --output-path {output.rep_seqs_legacy}
            rm -f {params.fasta}
            """

    rule Silva138_V3V4_taxonomy:
        input:
            rep_seqs = QIIME_DIR / "rep-seqs-analysis-legacy.qza",
            silva138_nb_classifier = REF_DIR / "Silva138" / "silva138_AB_V3-V4_classifier.qza"
        output:
            taxonomy = QIIME_DIR / "Silva138-taxonomy.qza",
            taxonomy_qzv = QIIME_DIR / "Silva138-taxonomy.qzv",
            sentinel = QIIME_DIR / ".Silva138_taxonomy_done"
        container:
            SILVA138_LEGACY_CLASSIFY_CONTAINER
        resources:
            mem_mb = 20000
        message:
            "[Silva138] Classifying 16S region V3V4 sequences using Silva138 Naive Bayes Classifier (legacy scikit-learn 0.24.1 container)..."
        shell:
            """
            set -euo pipefail
            qiime feature-classifier classify-sklearn \
                --i-classifier {input.silva138_nb_classifier} \
                --i-reads {input.rep_seqs} \
                --o-classification {output.taxonomy}
            echo "[Silva138] Generating taxonomy visualization..."
            qiime metadata tabulate \
                --m-input-file {output.taxonomy} \
                --o-visualization {output.taxonomy_qzv}
            touch {output.sentinel}
            """

    rule Silva138_collapse_to_genus:
        input:
            table    = QIIME_DIR / "table-analysis.qza",
            taxonomy = QIIME_DIR / "Silva138-taxonomy.qza",
            sentinel = QIIME_DIR / ".Silva138_taxonomy_done"
        output:
            genus_table = QIIME_DIR / "Silva138-genus-table-collapsed.qza",
            sentinel    = QIIME_DIR / ".Silva138_genus_collapsed_done"
        container:
            QIIME_CONTAINER
        resources:
            mem_mb = 4000
        message:
            "[Silva138] Collapsing NGS ASV table to genus level (Silva138 rank 6)..."
        shell:
            """
            qiime taxa collapse \
                --i-table {input.table} \
                --i-taxonomy {input.taxonomy} \
                --p-level 6 \
                --o-collapsed-table {output.genus_table}
            touch {output.sentinel}
            """

    rule Silva138_export_genus_biom:
        input:
            genus_table_qza = QIIME_DIR / "Silva138-genus-table-collapsed.qza",
            sentinel        = QIIME_DIR / ".Silva138_genus_collapsed_done"
        output:
            genus_biom = TABLES_DIR / "Silva138-study-seqs-genus.biom"
        container:
            QIIME_CONTAINER
        resources:
            mem_mb = 2000
        message:
            "[Silva138] Exporting NGS genus-collapsed feature table BIOM format..."
        shell:
            """
            qiime tools export \
                --input-path {input.genus_table_qza} \
                --output-path {TABLES_DIR}/exported_genus_temp_silva138
            mv {TABLES_DIR}/exported_genus_temp_silva138/feature-table.biom {output.genus_biom}
            rm -r {TABLES_DIR}/exported_genus_temp_silva138
            """

    rule Silva138_export_genus_taxonomy_as_tsv:
        input:
            taxonomy_qza = QIIME_DIR / "Silva138-taxonomy.qza",
            sentinel     = QIIME_DIR / ".Silva138_taxonomy_done"
        output:
            genus_taxonomy_tsv = TABLES_DIR / "exported-taxonomy" / "Silva138_genus_taxonomy.tsv"
        container:
            QIIME_CONTAINER
        resources:
            mem_mb = 2000
        message:
            "[Silva138] Exporting genus-collapsed taxonomy BIOM to TSV..."
        shell:
            """
            qiime tools export \
                --input-path {input.taxonomy_qza} \
                --output-path {TABLES_DIR}/exported-taxonomy-genus-temp-silva138
            mv {TABLES_DIR}/exported-taxonomy-genus-temp-silva138/taxonomy.tsv \
               {output.genus_taxonomy_tsv}
            rm -r {TABLES_DIR}/exported-taxonomy-genus-temp-silva138
            """

elif region == "full_length":
    # Trained via RESCRIPt against SILVA 138.2 (newer than the 138 release used by
    # QIIME2's own official classifier), with dereplicate + cull_seqs preprocessing
    # applied before fitting — RESCRIPt's recommended practice for classifier
    # accuracy. Built against scikit-learn 1.4.2, matching QIIME_CONTAINER, so no
    # legacy container is needed. Also replaces the lab's old
    # silva138_noEuk_AB_classifier.qza, which was trained on scikit-learn 0.24.1 and
    # can no longer be unpickled by the current environment.
    rule Silva138_Full_taxonomy:
        input:
            rep_seqs = QIIME_DIR / "rep-seqs-analysis.qza",
            silva138_nb_classifier = REF_DIR / "Silva138" / "SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza"
        output:
            taxonomy = QIIME_DIR / "Silva138-taxonomy.qza",
            taxonomy_qzv = QIIME_DIR / "Silva138-taxonomy.qzv",
            sentinel = QIIME_DIR / ".Silva138_taxonomy_done"
        container:
            QIIME_CONTAINER
        resources:
            mem_mb = 20000
        message:
            "[Silva138] Classifying 16S Full-length sequences using Silva138 Naive Bayes Classifier..."
        shell:
            """
            qiime feature-classifier classify-sklearn \
                --i-classifier {input.silva138_nb_classifier} \
                --i-reads {input.rep_seqs} \
                --o-classification {output.taxonomy}
            echo "[Silva138] Generating taxonomy visualization..."
            qiime metadata tabulate \
                --m-input-file {output.taxonomy} \
                --o-visualization {output.taxonomy_qzv}
            touch {output.sentinel}
            """

# Unlike Greengenes2 (which ships a static backbone tree and maps ASVs onto it via
# the greengenes2 plugin), Silva138 has no bundled backbone tree. Its phylogeny is
# built per-study by inserting ASVs into a SILVA reference alignment/tree with
# fragment-insertion (SEPP), then the table/rep-seqs are filtered down to only the
# features that were successfully placed. Here we use the Silva128 (not Silva138)
# reference tree/alignment because it's the only one available in QIIME2's SEPP
# plugin. View this post on the QIIME2 forum for more details:
# https://forum.qiime2.org/t/compatibility-of-sepp-silva128-with-taxonomy-classification-using-silva138/31172/3

rule Silva138_phylogeny:
    input:
        table = QIIME_DIR / "table-analysis.qza",
        rep_seqs = QIIME_DIR / "rep-seqs-analysis.qza",
        sepp_reference = REF_DIR / "Silva138" / "sepp-refs-silva-128.qza",
        sentinel = QIIME_DIR / ".Silva138_taxonomy_done"
    output:
        tree = QIIME_DIR / "Silva138-tree.qza",
        placements = QIIME_DIR / "Silva138-placements.qza",
        phylogeny_table = QIIME_DIR / "Silva138-table.qza",
        removed_table = QIIME_DIR / "Silva138-removed-table.qza",
        phylogeny_rep_seqs = QIIME_DIR / "Silva138-rep-seqs.qza",
        sentinel = QIIME_DIR / ".Silva138_phylogeny_done"
    # SEPP's placement step (pplacer) spawns one memory-heavy worker per thread;
    # at n_threads=14 this OOM-killed on a 30GB host placing full-length reads
    # against the SILVA 128 tree. Capped independently of n_threads rather than
    # lowering it globally, since other rules aren't memory-bound the same way.
    threads: 3
    container:
        QIIME_CONTAINER
    # Already came within a hair of OOM at n_threads=14 on this 30GB host (see
    # threads comment above) — budget it near the top of what's available.
    # runtime: long-running by design (see scripts/watch_sepp.sh, built
    # specifically to monitor this rule) -- had no explicit ceiling before,
    # silently inheriting whatever the SLURM profile's default-resources
    # runtime happened to be (previously 60min, which would have killed it).
    resources:
        mem_mb = 24000,
        runtime = 720
    message:
        "[Silva138] Building phylogenetic tree for 16S (V3V4 or Full-length) sequences via fragment-insertion (SEPP)..."
    shell:
        """
        set -euo pipefail
        # The Snakefile's global shell.prefix() disables Snakemake's default
        # set -euo pipefail for every rule, so a failure here (e.g. sepp being
        # OOM-killed) would otherwise silently fall through to filter-features,
        # which would then fail on a confusing "tree does not exist" error
        # instead of the real one. Set explicitly for this rule alone.

        # qiime's own --verbose gives nothing useful for this action: SEPP
        # partitions the reference tree into placement subsets (P_<n>), each
        # split into alignment subsets (A_<n>_<m>), and its underlying
        # run_sepp.py process logs each subset's progress and the totals it
        # computed up front — but that log is redirected to a file under /tmp
        # instead of reaching stdout/stderr, at a path randomly generated per
        # run (can't be hardcoded). Every 60s, find it by walking this shell's
        # process tree via /proc, then print one line: timestamp, job type,
        # which subset, and completed/total subsets — plus aggregate resource
        # usage across the whole process subtree. See scripts/sepp_monitor_lib.sh
        # (shared with scripts/watch_sepp.sh, which can attach to a run that's
        # already in flight from any terminal).
        source {WORKFLOW_DIR}/scripts/sepp_monitor_lib.sh
        (
            sepp_log=""
            placement_total="?"
            alignment_total="?"
            while true; do
                sleep 60
                if [ -z "$sepp_log" ]; then
                    sepp_log=$(find_sepp_log $$ || true)
                fi
                if [ -n "$sepp_log" ] && [ -r "$sepp_log" ]; then
                    if [ "$placement_total" = "?" ] || [ "$alignment_total" = "?" ]; then
                        read -r placement_total alignment_total <<< "$(sepp_totals "$sepp_log")"
                    fi
                    format_sepp_progress "$sepp_log" "$placement_total" "$alignment_total" >&2
                else
                    echo "[Silva138_phylogeny] still in setup — no SEPP log yet" >&2
                fi
                pids="$$ $(collect_descendants $$)"
                resource_summary "$pids" >&2
            done
        ) &
        monitor_pid=$!
        trap 'kill "$monitor_pid" 2>/dev/null' EXIT

        qiime fragment-insertion sepp \
            --i-representative-sequences {input.rep_seqs} \
            --i-reference-database {input.sepp_reference} \
            --p-threads {threads} \
            --p-alignment-subset-size 1000 \
            --p-placement-subset-size 5000 \
            --verbose \
            --o-tree {output.tree} \
            --o-placements {output.placements}
        qiime fragment-insertion filter-features \
            --i-table {input.table} \
            --i-tree {output.tree} \
            --o-filtered-table {output.phylogeny_table} \
            --o-removed-table {output.removed_table}
        qiime feature-table filter-seqs \
            --i-data {input.rep_seqs} \
            --i-table {output.phylogeny_table} \
            --o-filtered-data {output.phylogeny_rep_seqs}
        touch {output.sentinel}
        """

rule Silva138_taxa_barplot_qiime:
    input:
        table = QIIME_DIR / "table-analysis.qza",
        taxonomy = QIIME_DIR / "Silva138-taxonomy.qza",
        metadata = metadata_path
    output:
        barplot_qzv = QIIME_DIR / "Silva138-taxa-bar-plots.qzv",
        sentinel = QIIME_DIR / ".Silva138_taxa_barplot_done"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 4000
    message:
        "[Silva138] Generating QIIME taxa barplot visualization..."
    shell:
        """
        qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output.barplot_qzv}
        touch {output.sentinel}
        """

rule Silva138_export_taxonomy_as_tsv:
    input:
        taxonomy_qza = QIIME_DIR / "Silva138-taxonomy.qza"
    output:
        taxonomy = TABLES_DIR / "exported-taxonomy" / "Silva138_taxonomy.tsv"
    container:
        QIIME_CONTAINER
    resources:
        mem_mb = 2000
    message:
        "[Silva138] Exporting Silva138 taxonomy as TSV..."
    shell:
        """
        qiime tools export \
            --input-path {input.taxonomy_qza} \
            --output-path {TABLES_DIR}/exported-taxonomy-temp-silva138

        # Rename exported taxonomy file to your expected filename
        mv {TABLES_DIR}/exported-taxonomy-temp-silva138/taxonomy.tsv \
           {output.taxonomy}
        rm -r {TABLES_DIR}/exported-taxonomy-temp-silva138
        """
