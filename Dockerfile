FROM condaforge/miniforge3:24.3.0-0

WORKDIR /opt/snakemake_usearch/
COPY . .
WORKDIR /wd
ENV CONDA_DIR=/opt/conda

RUN export DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  # Remove imagemagick due to https://security-tracker.debian.org/tracker/CVE-2019-10131
  && apt-get purge -y imagemagick imagemagick-6-common \
  # Install common packages, non-root user
  # && apt-get install xxxx \
  && apt-get autoremove -y && apt-get clean -y && rm -rf /var/lib/apt/lists/*

RUN conda env create -f /opt/snakemake_usearch/environment.yml -n snakemake_usearch \
  && echo "conda activate snakemake_usearch" >> ~/.bashrc 

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "snakemake_usearch", "/bin/bash", "-c"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake_usearch"]
