FROM condaforge/miniforge3:24.3.0-0

WORKDIR /opt/snakemake_usearch/
COPY . .
WORKDIR /wd
# ideally also user a CONDA_ENV variable, but it doesn't expand in SHELL and ENTRYPOINT
ENV CONDA_DIR=/opt/conda

RUN export DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  # Remove imagemagick due to https://security-tracker.debian.org/tracker/CVE-2019-10131
  && apt-get purge -y imagemagick imagemagick-6-common \
  # Install common packages, non-root user
  # && apt-get install xxxx \
  && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# create conda environment
RUN conda env create -f /opt/snakemake_usearch/environment.yml -n snakemake_usearch \
  && echo "conda activate snakemake_usearch" >> ~/.bashrc 

# install usearch11 into the conda environment (or overwrite usearch12)
RUN wget https://github.com/rcedgar/usearch_old_binaries/raw/refs/heads/main/bin/usearch11.0.667_i86linux64 -O "/opt/conda/envs/snakemake_usearch/bin/usearch" \
  && chmod +x "/opt/conda/envs/snakemake_usearch/bin/usearch"

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "snakemake_usearch", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake_usearch"]
