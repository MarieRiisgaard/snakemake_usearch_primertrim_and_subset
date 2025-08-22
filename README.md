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
The pre-built Docker container available from [`ghcr.io/kasperskytte/snakemake_usearch:main`](https://github.com/KasperSkytte/snakemake_usearch/pkgs/container/snakemake_usearch) is built from the `Dockerfile` and includes all required software (including `usearch` version 11, not 12, see below).

### Conda
Requirements are listed in `environment.yml`. To create as a conda environment simply run:
```
conda env create --file environment.yml -n snakemake_usearch
```

All rules use the same environment, hence it's best to just load the environment first before running `snakemake`, and avoid using `--use-conda` or `--sdm conda`.

### usearch version 12 only includes a minimal subset of commands
The `usearch` version currently available (version 12) from the bioconda channel only contains a small subset of commands and doesn't include some required commands, for example `fastx_relabel` and `otutab_rare`. As `usearch` is a single binary requiring no external dependencies, you can install `usearch` version 11 by simply overwriting the `usearch` binary in the conda environment by downloading one of the precompiled binaries from [usearch_old_binaries](https://github.com/rcedgar/usearch_old_binaries/) matching your platform and architecture, and then manually place it in the `bin` folder in the conda environment (remember to make it executable). For example, on Linux, you can run the following:

```
conda env create -f environment.yml -n snakemake_usearch
conda activate snakemake_usearch
wget https://github.com/rcedgar/usearch_old_binaries/raw/refs/heads/main/bin/usearch11.0.667_i86linux64 -O "${CONDA_PREFIX}/bin/usearch"
chmod +x "${CONDA_PREFIX}/bin/usearch"
```

For other platforms or architectures adjust accordingly.

