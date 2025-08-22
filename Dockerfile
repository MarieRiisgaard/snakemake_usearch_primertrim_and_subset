FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2e184525363f6bdd065020f5415d2b4bc2a1cdc81dfa78829023880484d0ac0f"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/snakemake_usearch.yml
#   prefix: /conda-envs/389c17c98ac4659ba7ff04e3da4b7a48
#   name: snakemake_usearch
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - gzip=1.13
#     - filtlong=0.2.1
#     - cutadapt=5.1
#     - usearch=12.0_beta
RUN mkdir -p /conda-envs/389c17c98ac4659ba7ff04e3da4b7a48
COPY workflow/envs/snakemake_usearch.yml /conda-envs/389c17c98ac4659ba7ff04e3da4b7a48/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/389c17c98ac4659ba7ff04e3da4b7a48 --file /conda-envs/389c17c98ac4659ba7ff04e3da4b7a48/environment.yaml && \
    conda clean --all -y

# install usearch11 into the conda environment (or overwrite usearch12)
RUN wget https://github.com/rcedgar/usearch_old_binaries/raw/refs/heads/main/bin/usearch11.0.667_i86linux64 -O "/conda-envs/389c17c98ac4659ba7ff04e3da4b7a48/bin/usearch" \
  && chmod +x "/conda-envs/389c17c98ac4659ba7ff04e3da4b7a48/bin/usearch"
