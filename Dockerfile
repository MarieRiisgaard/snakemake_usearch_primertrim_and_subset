FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="fee386467746a52973fc78a834a9ea47ac2d46551f3f721de3594453d1c56380"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/snakemake_usearch.yml
#   prefix: /conda-envs/04a2b4c5b47798998d6f00ee3422dd77
#   name: snakemake_usearch
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - gzip=1.13
#     - filtlong=0.2.1
#     - cutadapt=5.1
#     - usearch=12.0_beta
RUN mkdir -p /conda-envs/04a2b4c5b47798998d6f00ee3422dd77
COPY workflow/envs/snakemake_usearch.yml /conda-envs/04a2b4c5b47798998d6f00ee3422dd77/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/04a2b4c5b47798998d6f00ee3422dd77 --file /conda-envs/04a2b4c5b47798998d6f00ee3422dd77/environment.yaml && \
    conda clean --all -y

# install usearch11 into the conda environment (or overwrite usearch12)
RUN wget https://github.com/rcedgar/usearch_old_binaries/raw/refs/heads/main/bin/usearch11.0.667_i86linux64 -O "/conda-envs/389c17c98ac4659ba7ff04e3da4b7a48/bin/usearch" \
  && chmod +x "/conda-envs/389c17c98ac4659ba7ff04e3da4b7a48/bin/usearch"
