rule concat_all:
    input:
        expand(os.path.join(config["tmp_dir"], "01-sample_prep", "{sample}", "{sample}_filtered_renamed.fastq"), sample=sample_dirs)
    output:
        temp(os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed.fastq"))
    log:
        os.path.join(config["log_dir"], "02-denoise", "concat_all.log")
    message:
        "Concatenating all samples before denoising (generate zOTUs/ASVs)",
    resources:
        mem_mb=512,
        runtime=30,
        cpus_per_task=1
    threads: 1
    shell:
      """
        exec &> "{log}"
        set -euxo pipefail
        
        cat {input} > {output}
      """

rule orient:
    input:
        os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed.fastq")
    output:
        fq=temp(os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed_oriented.fastq")),
        tab=temp(os.path.join(config["tmp_dir"], "02-denoise", "orient.txt"))
    log:
        os.path.join(config["log_dir"], "02-denoise", "orient.log")
    params:
        db=config["db_sintax"]
    message:
        "Orienting reads"
    resources:
        mem_mb=2048,
        runtime=120,
        cpus_per_task=1
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail

        usearch11 -orient {input} \
          -db {params.db} \
          -fastqout {output.fq} \
          -tabbedout {output.tab}
        """

rule trim_reads:
    input:
        fq=os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed_oriented.fastq")
    output:
        fq=os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed_oriented_trimmed.fastq")
    log:
        os.path.join(config["log_dir"], "02-denoise", "trim_reads.log")
    params:
        truncate_length=config["truncate_length"]
    message:
        "Trimming oriented reads"
    resources:
        mem_mb=1024,
        runtime=60,
        cpus_per_task=1
    threads: 1
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail

        usearch -fastx_truncate {input.fq} \
          -trunclen {params.truncate_length} \
          -fastqout {output.fq}
        """

# according to Edgar QC filtering should be done here, not in sample_prep
# rule qc_filter:
#   shell:
#     """
#     usearch -fastq_filter trimmed.fq -fastq_maxee 1.0 -fastaout filtered.fa
#     """

rule derep:
    input:
      os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed_oriented_trimmed.fastq")
    output:
      temp(os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed_oriented_trimmed_derep.fa"))
    log:
      os.path.join(config["log_dir"], "02-denoise", "derep.log")
    message:
        "Dereplicating all samples"
    resources:
        mem_mb=lambda wc, input: max(1.5 * input.size_mb, 4096), # a 64GB file took 72GB mem
        runtime=600,
        cpus_per_task=4
    params:
      derep_minsize=config["derep_minsize"]
    threads: 4 # this command spends most of the time just reading in the file
    shell:
      """
        exec &> "{log}"
        set -euxo pipefail
        
        usearch -fastx_uniques \
          {input} \
          -fastaout {output} \
          -sizeout \
          -threads {threads} \
          -minuniquesize {params.derep_minsize} \
          -relabel Uniq
      """

rule unoise:
    input:
      os.path.join(config["tmp_dir"], "02-denoise", "all_samples_filtered_renamed_oriented_trimmed_derep.fa")
    output:
      os.path.join(config["output_dir"], "zOTUs.fa")
    log:
      os.path.join(config["log_dir"], "02-denoise", "unoise.log")
    message:
      "Denoising/generating zOTUs/ASVs"
    resources:
      mem_mb=lambda wc, input: max(3 * input.size_mb, 512),
      runtime=600,
      cpus_per_task=1
    params:
      unoise_minsize=config["unoise_minsize"]
    threads: 1
    shell:
      """
        exec &> "{log}"
        set -euxo pipefail
        
        usearch -unoise3 \
          {input} \
          -zotus {output} \
          -sizeout \
          -minsize {params.unoise_minsize}
    """
