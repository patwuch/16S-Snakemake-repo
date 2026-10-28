# scripts/build_manifest.py
import os
# Send this job's stdout/stderr, subprocesses included, to its Snakemake log.
_log = os.open(snakemake.log[0], os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
os.dup2(_log, 1)
os.dup2(_log, 2)

import glob
import pandas as pd

# Extract Snakemake inputs as strings
raw_dir = snakemake.params.raw_dir
sequence_type = snakemake.params.sequence_type
region = snakemake.params.region
study_dir = snakemake.params.study_dir
manifest_file = snakemake.output.manifest_file

manifest_data = []

# Loop through immediate subdirectories of raw_dir
for sample_id in os.listdir(raw_dir):
    sample_folder = os.path.join(raw_dir, sample_id)
    if not os.path.isdir(sample_folder):
        continue

    files = glob.glob(os.path.join(sample_folder, "*"))
    print(f"Files found for {sample_id}: {[os.path.basename(f) for f in files]}")

    if sequence_type == "NGS" and region == "region_V3V4":  # paired-end
        forward_file = None
        reverse_file = None

        for f in files:
            name = os.path.basename(f)
            if name.endswith("_1.fq.gz") or name.endswith("_1.fastq.gz") or name.endswith("_R1.fastq.gz"):
                forward_file = os.path.abspath(f)
            elif name.endswith("_2.fq.gz") or name.endswith("_2.fastq.gz") or name.endswith("_R2.fastq.gz"):
                reverse_file = os.path.abspath(f)

        if not forward_file or not reverse_file:
            raise ValueError(f"Missing forward or reverse read for sample {sample_id}")

        manifest_data.append({
            "sample-id": sample_id,
            "forward-absolute-filepath": forward_file,
            "reverse-absolute-filepath": reverse_file
        })

    elif sequence_type == "TGS" and region == "full_length":  # single-end
        for f in files:
            manifest_data.append({
                "sample-id": sample_id,
                "absolute-filepath": os.path.abspath(f)
            })

    else:
        raise ValueError("Unsupported sequence type or region configuration.")

# Write manifest
manifest_path = os.path.join(study_dir, "manifest.tsv")
pd.DataFrame(manifest_data).to_csv(manifest_path, index=False, sep='\t')
print(f"Manifest written to {manifest_path}")
