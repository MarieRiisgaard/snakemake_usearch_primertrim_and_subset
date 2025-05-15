rule concat_all:
    input:
        expand(os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}_filtered_renamed.fastq"), sample=sample_dirs)
    output:
        temp(os.path.join(config["tmp_dir"], "samples", "all_samples_filtered_renamed.fastq"))
    log:
        os.path.join(config["log_dir"], "concat_all.log")
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

#rule orient:

rule derep:
    input:
      os.path.join(config["tmp_dir"], "samples", "all_samples_filtered_renamed.fastq")
    output:
      temp(os.path.join(config["tmp_dir"], "samples", "all_samples_filtered_renamed_derep.fa"))
    log:
      os.path.join(config["log_dir"], "derep.log")
    conda:
        "../envs/env.yml"
    message:
        "Dereplicating all samples before denoising (generate zOTUs/ASVs)"
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
      os.path.join(config["tmp_dir"], "samples", "all_samples_filtered_renamed_derep.fa")
    output:
      os.path.join(config["output_dir"], "zOTUs.fa")
    log:
      os.path.join(config["log_dir"], "unoise.log")
    conda:
      "../envs/env.yml"
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
