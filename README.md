# Core Laboratory of Human Microbiome @ Taipei Medical University

Our lab provides end-to-end services from storing biological samples to publication-ready reports. This includes step like DNA extraction, artefact creation, taxa classification, diversity + phylogenetic analysis, and functional metabolomics prediction. Currently we handle:

* Next Generation Sequencing (16S V3-V4 rRNA)
* Third Generation Sequencing (16S Full Length)
* Shotgun Metagenomics
* Anaerobic Bacterium Cultivation

For more details, please visit [our landing page](https://microbiome-in-tmu.mystrikingly.com).

## Table of Contents

* __main__: Store your project specific data here by folder
* __reference__: Store your taxonomy or phylogenetic reference db here
* __config__: Store your project specific analysis configuration file (.yaml) here. Match name with project folder in __main__

See [`quickstart.md`](quickstart.md) for metadata.tsv column-naming rules and general Snakemake usage tips (dry runs, `--until`, `--touch`, etc), and [`pipeline_overview.md`](pipeline_overview.md) for a stage-by-stage walkthrough of what the Snakemake workflow actually does, from manifest building through the final PDF report.

## Currently Supported Taxonomic/Phylogeny Databases

* Greengenes2
* SILVA 138 (see [`pipeline_overview.md`](pipeline_overview.md) Stage 3b for its SEPP phylogeny caveats)

## Currently Supported Functional Prediction Database

* PICRUST2

## Currently Supported Differential Analysis Methods

* LEfSe (the original Huttenhower implementation — see [`pipeline_overview.md`](pipeline_overview.md) Stage 5)

Some projects also use other differential abundance methods, like MaAsLin3, ALDEx3 and ANCOM-BC2. Those aren't part of the standard run; you'll find what each one does in that project's own folder.

## Naming Convention Guide

It is suggested to keep your naming convention simple, i.e. use hyphens for QIIME2 artefacts(.qza, .qzv) and underscores for everything else (.tsv, .csv) so you know which files are meant to be only processed via QIIME2.

## Environment Guide

This is the Linux repository. Per-rule software (QIIME2, PICRUSt2, LEfSe, the report compiler) runs inside Apptainer containers rather than ad hoc conda environments, so a run is reproducible regardless of what happens to be installed on the host.

**Running the pipeline**: activate an environment with Snakemake and Apptainer installed (`environment.yml` defines this — e.g. `conda env create -f environment.yml && conda activate smk9`), then invoke Snakemake with `--sdm apptainer`:

```Shell
snakemake --sdm apptainer --configfile config/<your_config>.yaml <target>
```

`--sdm` (`--software-deployment-method`) tells Snakemake to run each rule's `container:` directive through Apptainer instead of provisioning a conda env per rule.

**Prefer running through `scripts/run_snakemake.sh` instead of calling `snakemake` directly** — it's a drop-in wrapper (same arguments) that runs the pipeline inside a transient systemd cgroup with swap disabled, so a rule that runs away on memory gets killed fast instead of thrashing the whole host into swap:

```Shell
scripts/run_snakemake.sh --sdm apptainer --configfile config/<your_config>.yaml <target>
```

Override the memory ceiling with `SNAKEMAKE_MEM_MAX` (cgroup cap, default `26G`) and `SNAKEMAKE_MEM_BUDGET` (the `mem_mb` resource budget passed to Snakemake, default `24000`).

Before starting, the wrapper runs `scripts/check_containers.sh`, which refuses to run if any `.sif` the workflow uses no longer matches its `.def` or `envs/*.yaml`. It doesn't rebuild anything; rebuild with `containers/build.sh --force <name>`, or set `SKIP_CONTAINER_CHECK=1` to run anyway.

**Building/rebuilding containers**: each `containers/<name>.sif` is built from the matching `containers/<name>.def`, which in turn installs the matching `envs/<name>.yaml` conda spec. So in plain English: `envs/*.yaml` stays the source of truth for package lists, want to change something? Edit those, not a built `.sif` directly. Just run `containers/build.sh` to build any `.def` missing a `.sif`, or `containers/build.sh --force <name>` to rebuild one after editing its `envs/*.yaml`. Since we don't use conda anymore, snakemake does not rebuild containers on its own, and builds must run one at a time (`build.sh` already does this) — concurrent `apptainer build --fakeroot` runs corrupt each other via fakeroot namespace contention.

**Filesystem access inside containers**: Snakemake invokes Apptainer without `--contain`/`--containall`, so containers share the host filesystem by default — `main/` and `reference/` are already visible inside every container with no `--bind` configuration needed. Do not add `--contain`/`--containall` via `--apptainer-args` unless you also add explicit `--bind` flags for `main/` and `reference/`, or those paths will silently disappear inside the container.

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
