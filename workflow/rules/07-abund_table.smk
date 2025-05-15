# MUUUCH faster to do this per sample and merge afterwards
# usearch -otutab does NOT scale linearly with more threads
rule abund_table:
    input:
        zotus=os.path.join(config["output_dir"], "zOTUs.fa"),
        allreads_unfiltered=os.path.join(config["tmp_dir"], "samples", "{sample}", "{sample}_renamed.fastq")
    output:
        temp(os.path.join(config["tmp_dir"], "abund_table", "{sample}_abund_table.tsv"))
    log:
        os.path.join(config["log_dir"], "abund_table", "{sample}_abund_table.log")
    message:
        "Estimating abundances of zOTUs/ASVs in each sample",
    resources:
        mem_mb=512,
        runtime=120,
        cpus_per_task=4
    threads: 4
    params:
        sample_sep=config['sample_sep']
    conda:
        "../envs/env.yml"
    shell:
        """
            exec &> "{log}"
            set -euxo pipefail
        
            usearch -otutab \
                "{input.allreads_unfiltered}" \
                -zotus "{input.zotus}" \
                -otutabout "{output}" \
                -threads "{threads}" \
                -sample_delim "{params.sample_sep}"
        """

rule merge_abund_tables:
    input:
        expand(os.path.join(config["tmp_dir"], "abund_table", "{sample}_abund_table.tsv"), sample=sample_dirs)
    output:
        os.path.join(config["output_dir"], "abund_table.tsv")
    log:
        os.path.join(config["log_dir"], "merge_abund_tables.log")
    message:
        "Merging abundance tables",
    resources:
        mem_mb=2048,
        runtime=60,
        cpus_per_task=1
    threads: 1
    conda:
        "../envs/env.yml"
    params:
        input_csv=lambda wildcards, input: ",".join(input)
    shell:
        """
            exec &> "{log}"
            set -euxo pipefail
        
            # merge abundance tables
            usearch11 -otutab_merge \
                {params.input_csv} \
                -output "{output}"
        """