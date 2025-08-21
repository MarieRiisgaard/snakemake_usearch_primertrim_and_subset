# snakemake_usearch
Simple snakemake pipeline based on the standard [usearch11 SOP](https://drive5.com/usearch/manual/uparse_pipeline.html).

## Usage
First install snakemake and the required software into a conda environment or use the container as described below. Then run using fx:
```
conda activate snakemake_usearch
snakemake --jobs 96
```

Use an executor if you are running on a HPC cluster. See the `slurm_submit.sbatch` for an example when running on a SLURM cluster.

## Requirements
Install the required software by using either the provided `Dockerfile` or `environment.yml` file to build a Docker container or conda environment with all the required tools, see below.

### Docker
Tre-built Docker container available from [`ghcr.io/kasperskytte/snakemake_usearch:main`](https://github.com/KasperSkytte/snakemake_usearch/pkgs/container/snakemake_usearch) includes all required software, including `usearch` version 11.

### Conda
Requirements are listed in `environment.yml`. To create as a conda environment simply run:
```
conda env create --file environment.yml -n snakemake_usearch
```

#### usearch version 12 is insiffucient
The `usearch` version currently available (version 12) from the bioconda channel only contains a small subset of commands and doesn't include `-otutab_rare`, which is required. If  you don't need a rarefied abundance table you can remove the `rarefy_abund_table` rule and version 12 would be sufficient. Otherwise you must download a precompiled `usearch` version 11 binary from [usearch_old_binaries](https://github.com/rcedgar/usearch_old_binaries/) and manually place it in the `bin` folder in the conda environment and make it executable. For example, on Linux, you can run the following:

```
conda env create -f environment.yml -n snakemake_usearch
conda activate snakemake_usearch
wget https://github.com/rcedgar/usearch_old_binaries/raw/refs/heads/main/bin/usearch11.0.667_i86linux64 -O "${CONDA_PREFIX}/bin/usearch"
chmod +x "${CONDA_PREFIX}/bin/usearch"
```

For other platforms or architectures adjust accordingly.

All rules use the same environment, hence don't start the workflow with `--use-conda`, instead just load the environment first before running `snakemake`.
