############################################
# 05 – Generating of SINTAX taxonomy predictions
############################################
rule sintax_subset:
    input:
        os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa")
    output:
        os.path.join(config["output_dir"], "05-sintax", "{subset}", "zOTUs.sintax")
    log:
        os.path.join(config["log_dir"], "05-sintax", "{subset}_sintax.log")
    message:
        "Predicting taxonomy of zOTUs for subset {wildcards.subset} using SINTAX"
    resources:
        mem_mb=4096,
        runtime=60,
        cpus_per_task=config["max_threads"]
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: config["max_threads"]
    params:
        db=config["db_sintax"]
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p $(dirname {output})

        usearch -sintax \
          "{input}" \
          -db "{params.db}" \
          -tabbedout "{output}" \
          -strand both \
          -sintax_cutoff 0.8 \
          -threads {threads}

        sort -V "{output}" -o "{output}"

        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
        """
