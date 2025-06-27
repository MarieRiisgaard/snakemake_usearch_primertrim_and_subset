# snakemake_usearch

## Installing required software
Requirements are listed in `environment.yml`. To create as a conda environment simply run:
```
conda env create --file environment.yml -n snakemake_usearch
```

and then install usearch11 manually into the environment afterwards as described below.

### Installing usearch11 into the conda environment
The open source [`usearch12`](https://github.com/rcedgar/usearch12) does not contain all the required commands, and `usearch11` is currently not available from any conda channels, so you must download a precompiled usearch11 binary from [usearch_old_binaries](https://github.com/rcedgar/usearch_old_binaries/) and manually place it in the `bin` folder in the conda environment. For example, on Linux, you can run the following:

```
conda env create -f environment.yml -n snakemake_usearch
conda activate snakemake_usearch
wget https://github.com/rcedgar/usearch_old_binaries/raw/refs/heads/main/bin/usearch11.0.667_i86linux64 -O "${CONDA_PREFIX}/bin/usearch"
chmod +x "${CONDA_PREFIX}/bin/usearch"
```

For other platforms or architectures adjust accordingly. All rules use the same environment, hence don't start the workflow with `--use-conda`.
