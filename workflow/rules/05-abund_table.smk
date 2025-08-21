# MUUUCH faster to do this per sample and merge afterwards
# usearch -otutab does NOT scale linearly with more threads
rule abund_table:
    input:
        zotus=os.path.join(config["output_dir"], "zOTUs.fa"),
        allreads_unfiltered=os.path.join(config["tmp_dir"], "01-sample_prep", "{sample}", "{sample}_renamed.fastq")
    output:
        temp(os.path.join(config["tmp_dir"], "abund_table", "abund_table_{sample}.tsv"))
    log:
        os.path.join(config["log_dir"], "05-abund_table", "abund_table_{sample}.log")
    message:
        "Estimating abundances of zOTUs/ASVs in each sample",
    resources:
        mem_mb=1024, # this needs to be calculated dynamically based on input file size
        runtime=120,
        cpus_per_task=4
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    threads: 4
    params:
        sample_sep=config['sample_sep']
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
            # sometimes output files are empty and this errors, what to do?!
        """

rule merge_abund_tables:
    input:
        expand(os.path.join(config["tmp_dir"], "abund_table", "abund_table_{sample}.tsv"), sample=sample_dirs)
    output:
        os.path.join(config["output_dir"], "abund_table.tsv")
    log:
        os.path.join(config["log_dir"], "05-abund_table", "merge_abund_tables.log")
    message:
        "Merging abundance tables",
    resources:
        mem_mb=2048,
        runtime=60,
        cpus_per_task=1
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    threads: 1
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
            # this errors if just one file is empty, what to do?!
        """

rule rarefy_abund_table:
    input:
        os.path.join(config["output_dir"], "abund_table.tsv")
    output:
        os.path.join(config["output_dir"], "abund_table_rarefied.tsv")
    log:
        os.path.join(config["log_dir"], "05-abund_table", "rarefy_abund_tables.log")
    message:
        "Rarefying abundance table"
    resources:
        mem_mb=2048,
        runtime=120,
        cpus_per_task=1
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    threads: 1
    params:
        rarefy=config['rarefy']
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail

        usearch -otutab_rare {input} -sample_size {params.rarefy} -output {output}
        """
