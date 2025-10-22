# List of all subsets (include all_reads + sizes)
subset_dirs = ["all_reads"] + [f"sample_size_{s}" for s in config.get("subsample_sizes", [10000, 20000, 50000, 100000, 200000, 400000])]

############################################
# 1️⃣ Concatenate all barcodes per subset
############################################

# --- Helper function to select correct inputs ---
def inputs_for_subset(wc):
    if wc.subset == "all_reads":
        # Read trimmed FASTQs directly from tmp/02-denoise
        return expand(
            os.path.join(
                config["tmp_dir"],
                "02-denoise",
                "{sample}",
                "{sample}_trimmed.fastq"
            ),
            sample=sample_dirs,
        )
    else:
        # Subsampled reads for other subsets
        size = wc.subset.replace("sample_size_", "")
        return expand(
            os.path.join(
                config["output_dir"],
                "subsample",
                wc.subset,
                "{sample}_subsampled_" + size + ".fastq"
            ),
            sample=sample_dirs,
        )


rule concat_subsample:
    input:
        inputs_for_subset
    output:
        temp(os.path.join(config["tmp_dir"], "04-denoise", "{subset}", "all_samples_trimmed.fastq"))
    log:
        os.path.join(config["log_dir"], "04-denoise", "{subset}_concat.log")
    message:
        "Concatenating all barcodes for subset {wildcards.subset}"
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 1
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p "$(dirname {output})"
        cat {input} > "{output}"

        # safety check
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
        """



############################################
# 2️⃣ Dereplicate per subset
############################################
rule derep_subset:
    input:
        os.path.join(config["tmp_dir"], "04-denoise", "{subset}", "all_samples_trimmed.fastq")
    output:
        temp(
            os.path.join(
                config["tmp_dir"],
                "04-denoise",
                "{subset}",
                "all_samples_trimmed_derep.fa"
            )
        )
    log:
        os.path.join(config["log_dir"], "04-denoise", "{subset}_derep.log")
    message:
        "Dereplicating subset {wildcards.subset}"
    params:
        derep_minsize=config["derep_minsize"]
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 4
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        usearch -fastx_uniques \
          {input} \
          -fastaout {output} \
          -sizeout \
          -threads {threads} \
          -minuniquesize {params.derep_minsize} \
          -relabel Uniq

        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
        """


############################################
# 3️⃣ UNOISE denoise per subset
############################################
rule unoise_subset:
    input:
        os.path.join(config["tmp_dir"], "04-denoise", "{subset}", "all_samples_trimmed_derep.fa")
    output:
        os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa")
    log:
        os.path.join(config["log_dir"], "04-denoise", "{subset}_unoise.log")
    message:
        "Denoising subset {wildcards.subset} with UNOISE3"
    params:
        unoise_minsize=config["unoise_minsize"]
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 1
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p $(dirname {output})

        usearch -unoise3 \
          {input} \
          -zotus {output} \
          -sizeout \
          -minsize {params.unoise_minsize}

        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
        """
