## Column Naming

For metadata.tsv formatting, first column is always 'sampleid'--no space, no hyphen, no capitalisation--just these exact lower case characters.

The second column, no matter what it is named is, automatically considered the primary grouping axis. I do advise naming it 'Group' for easier reference but it can really be anything.

The primary grouping axis helps you enforce both the order and the color of how group names appear (up-down, left-right) in plots.

By optionally creating a 'PrimaryOrder' and/or a 'PrimaryColor' column, plots generated will follow the integer (1,2,3...) and the color (#4D4D4D, #377EB8, #FFFFFF), so long as they don't violate the one-to-one matching between the primary grouping axis and specified values.

Non-primary grouping axis will be generated alphanumerically, and only one primary grouping axis can exist in a single run.

## Customising your config_project.yaml

Consult config_project.yaml in the config folder and follow its instructions.

## Starting a Snakemake Run

(A) Check you are in the right conda environment, run:

```Shell
conda activate smk9
```

(B) Make sure you are at the repo root, run

```Shell
cd /srv/projects/16S-Snakemake-repo
```

(C) Decide your snakemake run parameters and format them like below

```
--cores all --sdm apptainer --configfile config/config_TEMPLATE.yaml
```

(D) Append the string in step C to the helper script, so you run:

```Shell
scripts/run_snakemake.sh --cores all --sdm apptainer --configfile config/config_template.yaml
```

## Run as SLURM Batch Jobs

To launch multiple project runs in parallel, use **sbatch**.

```Shell
sbatch --job-name=tisha260707 scripts/slurm/submit_project.sbatch config/config_tisha260707.yaml
sbatch --job-name=hu260708 scripts/slurm/submit_project.sbatch config/config_hu260708.yaml
sbatch --job-name=ahmad260220 scripts/slurm/submit_project.sbatch config/config_ahmad260220
```

You can give --job-name any value you want, it does not necessarily need to match the config or the project folder as strictly. But for your own sake it is recommended to name it the same.

## Additional Tips

If you are worried whatever you run might accidentally trigger unwanted steps,
for safety add -n as one parameter to the string to perform a dry run so Snakemake will tell you what it wants to do. If you approve of all the steps, remove the -n and it will run exactly that.

Use --until to only run up to a certain rule

Use --touch to tell Snakemake to not rerun some thing but "pretend" it has been reran

Use --omit-from to not run a specific rule
