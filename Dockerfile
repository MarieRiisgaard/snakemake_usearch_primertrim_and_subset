FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="931e4f24e8b2de6086f8a79523a51a3c2825f268d1a4283b15406894653a5746"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/snakemake_usearch.yml
#   prefix: /conda-envs/b65a0fdd5bf9a7ad3ad6ede2fe067cc3
#   name: snakemake_usearch
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - gzip=1.13
#     - filtlong=0.2.1
#     - cutadapt=5.1
#     - usearch=11.0.667
RUN mkdir -p /conda-envs/b65a0fdd5bf9a7ad3ad6ede2fe067cc3
COPY workflow/envs/snakemake_usearch.yml /conda-envs/b65a0fdd5bf9a7ad3ad6ede2fe067cc3/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/b65a0fdd5bf9a7ad3ad6ede2fe067cc3 --file /conda-envs/b65a0fdd5bf9a7ad3ad6ede2fe067cc3/environment.yaml && \
    conda clean --all -y
