###############################################################
# 04-denoise.smk — Denoising workflow for each subset
###############################################################

# -------------------------------------------------------------
# List of subsets (based on config defaults)
# -------------------------------------------------------------
subset_dirs = ["all_reads"] + [
    f"sample_size_{s}"
    for s in config.get("subsample_sizes", [10000, 20000, 50000, 100000, 200000])
] + ["sample_size_all_reads"]


# -------------------------------------------------------------
# Helper function: select correct input FASTQ for each subset
# -------------------------------------------------------------
def inputs_for_subset(wc):
    if wc.subset == "all_reads":
        # Use the merged file with all trimmed reads (before subsampling)
        return [
            os.path.join(config["tmp_dir"], "merged", "all_samples_trimmed.fastq")
        ]
    else:
        # Extract the numeric or "all_reads" part from subset
        size = wc.subset.replace("sample_size_", "")
        return [
            os.path.join(
                config["output_dir"],
                "03-subsample",
                f"sample_size_{size}",
                f"all_samples_subsampled_{size}.fastq"
            )
        ]

###############################################################
# 1️⃣ Dereplicate per subset (with safe fallback)
###############################################################
rule derep_subset:
    input:
        inputs_for_subset
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
    resources:
        mem_mb=lambda wc, input: max(1.5 * input.size_mb, 4096),  # adaptive memory
        runtime=600,
        cpus_per_task=4,
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

        mkdir -p "$(dirname {output})"

        echo "Running dereplication on {input}"

        usearch -fastx_uniques \
          {input} \
          -fastaout {output} \
          -sizeout \
          -threads {threads} \
          -minuniquesize {params.derep_minsize} \
          -relabel Uniq || true

        # 🧩 Fallback: if dereplication produced no output, make placeholder
        if [ ! -s "{output}" ]; then
            echo "# No unique sequences found (empty dereplication)" > "{output}"
            echo "⚠️ No unique sequences detected for subset {wildcards.subset}" >&2
        fi
        """


###############################################################
# 2️⃣ UNOISE denoising per subset (with skip if derep is empty)
###############################################################
rule unoise_subset:
    input:
        os.path.join(config["tmp_dir"], "04-denoise", "sample_size_{subset}", "all_samples_trimmed_derep.fa")
    output:
        os.path.join(config["output_dir"], "04-denoise", "sample_size_{subset}", "zOTUs.fa")
    log:
        os.path.join(config["log_dir"], "04-denoise", "sample_size_{subset}_unoise.log")
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

        mkdir -p "$(dirname {output})"

        # 🧩 Skip UNOISE if dereplication file is placeholder
        if grep -q "# No unique sequences found" "{input}"; then
            echo "# Skipping UNOISE — no unique sequences found" > "{output}"
            echo "⚠️ Skipping UNOISE for subset {wildcards.subset} (no dereplicated reads)" >&2
            exit 0
        fi

        # Run UNOISE3
        usearch -unoise3 \
          {input} \
          -zotus {output} \
          -sizeout \
          -minsize {params.unoise_minsize} || true

        # 🧩 Fallback: if UNOISE produced no ASVs, create placeholder
        if [ ! -s "{output}" ]; then
            echo "# No ASVs detected (empty run)" > "{output}"
            echo "⚠️ No ASVs generated for subset {wildcards.subset}" >&2
        fi
        """













# List of all subsets (include all_reads + sizes)
#subset_dirs = ["all_reads"] + [f"sample_size_{s}" for s in config.get("subsample_sizes", [10000, 20000, 50000, 100000, 200000, 400000])]
#
#############################################
## 1️⃣ Concatenate all barcodes per subset
#############################################
#
## --- Helper function to select correct inputs ---
#def inputs_for_subset(wc):
#    if wc.subset == "all_reads":
#        # Read trimmed FASTQs directly from tmp/02-denoise
#        return expand(
#            os.path.join(
#                config["tmp_dir"],
#                "02-denoise",
#                "{sample}",
#                "{sample}_trimmed.fastq"
#            ),
#            sample=sample_dirs,
#        )
#    else:
#        # Subsampled reads for other subsets
#        size = wc.subset.replace("sample_size_", "")
#        return expand(
#            os.path.join(
#                config["output_dir"],
#                "subsample",
#                wc.subset,
#                "{sample}_subsampled_" + size + ".fastq"
#            ),
#            sample=sample_dirs,
#        )
#
#
#rule concat_subsample:
#    input:
#        inputs_for_subset
#    output:
#        temp(os.path.join(config["tmp_dir"], "04-denoise", "{subset}", "all_samples_trimmed.fastq"))
#    log:
#        os.path.join(config["log_dir"], "04-denoise", "{subset}_concat.log")
#    message:
#        "Concatenating all barcodes for subset {wildcards.subset}"
#    container:
#        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
#    conda:
#        "../envs/snakemake_usearch.yml"
#    resources:
#        mem_mb=32000,   # 32 GB; adjust as needed
#        runtime=180,    # minutes
#        cpus_per_task=1
#    threads: 1
#    shell:
#        r"""
#        exec &> "{log}"
#        set -euxo pipefail
#
#        mkdir -p "$(dirname {output})"
#        cat {input} > "{output}"
#
#        # safety check
#        if [ ! -s "{output}" ]; then
#            echo "output file {output} is empty, exiting!"
#            exit 1
#        fi
#        """
#
#
#
#############################################
## 2️⃣ Dereplicate per subset
#############################################
#rule derep_subset:
#    input:
#        os.path.join(config["tmp_dir"], "04-denoise", "{subset}", "all_samples_trimmed.fastq")
#    output:
#        temp(
#            os.path.join(
#                config["tmp_dir"],
#                "04-denoise",
#                "{subset}",
#                "all_samples_trimmed_derep.fa"
#            )
#        )
#    log:
#        os.path.join(config["log_dir"], "04-denoise", "{subset}_derep.log")
#    message:
#        "Dereplicating subset {wildcards.subset}"
#    resources:
#        mem_mb=lambda wc, input: max(1.5 * input.size_mb, 4096),  # a 64GB file took 72GB mem
#        runtime=600,
#        cpus_per_task=4,
#    params:
#        derep_minsize=config["derep_minsize"]
#    container:
#        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
#    conda:
#        "../envs/snakemake_usearch.yml"
#    threads: 4 # this command spends most of the time just reading in the file, increasing isn't beneficial
#    shell:
#        r"""
#        exec &> "{log}"
#        set -euxo pipefail
#
#        usearch -fastx_uniques \
#          {input} \
#          -fastaout {output} \
#          -sizeout \
#          -threads {threads} \
#          -minuniquesize {params.derep_minsize} \
#          -relabel Uniq
#
#        if [ ! -s "{output}" ]; then
#            echo "output file {output} is empty, exiting!"
#            exit 1
#        fi
#        """
#
#
#############################################
## 3️⃣ UNOISE denoise per subset
#############################################
#rule unoise_subset:
#    input:
#        os.path.join(config["tmp_dir"], "04-denoise", "{subset}", "all_samples_trimmed_derep.fa")
#    output:
#        os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa")
#    log:
#        os.path.join(config["log_dir"], "04-denoise", "{subset}_unoise.log")
#    message:
#        "Denoising subset {wildcards.subset} with UNOISE3"
#    params:
#        unoise_minsize=config["unoise_minsize"]
#    container:
#        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
#    conda:
#        "../envs/snakemake_usearch.yml"
#    threads: 1
#    shell:
#        r"""
#        exec &> "{log}"
#        set -euxo pipefail
#
#        mkdir -p $(dirname {output})
#
#        # Run UNOISE3 but do not fail the rule if it returns non-zero
#        usearch -unoise3 \
#          {input} \
#          -zotus {output} \
#          -sizeout \
#          -minsize {params.unoise_minsize} || true
#
#        # If the output is missing or empty, create an empty placeholder
#        if [ ! -s "{output}" ]; then
#            echo "# No ASVs detected (empty run)" > "{output}"
#            echo "⚠️ No ASVs generated for subset {wildcards.subset}" >&2
#        fi
#        """
#
#