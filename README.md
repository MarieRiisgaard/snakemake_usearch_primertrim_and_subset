# snakemake_usearch
Simple snakemake workflow based on the standard [usearch11 SOP](https://drive5.com/usearch/manual/uparse_pipeline.html).

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
The pre-built Docker container available from [`ghcr.io/kasperskytte/snakemake_usearch:main`](https://github.com/KasperSkytte/snakemake_usearch/pkgs/container/snakemake_usearch) is built from the `Dockerfile` and includes all required software.

### Conda
Requirements are listed in `environment.yml`. To create as a conda environment simply run:
```
conda env create --file environment.yml -n snakemake_usearch
```
