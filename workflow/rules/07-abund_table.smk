# MUUUCH faster to do this per sample and merge afterwards
# usearch -otutab does NOT scale linearly with more threads
# more than 16 has diminishing returns
rule abund_table:
    input:
        zotus=lambda wc: os.path.join(
            config["output_dir"], "04-denoise", f"sample_size_{wc.subset}", "zOTUs.fa.filtered"
        ),
        allreads_unfiltered=lambda wc: os.path.join(
            config["output_dir"], "03-subsample", f"sample_size_{wc.subset}", f"all_samples_subsampled_{wc.subset}.fastq"
        )
    output:
        os.path.join(config["tmp_dir"], "07-abund_table", "abund_table_{subset}.tsv"),
    log:
        os.path.join(config["log_dir"], "07-abund_table", "abund_table_{subset}.log"),
    message:
        "{wildcards.subset}: Estimating abundances of zOTUs/ASVs"
    threads: lambda wc, input: min(config.get("max_threads", 8), 16)
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    params:
        sample_sep=config.get("sample_sep", "_"),
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        if [ ! -s "{input.zotus}" ] || [ ! -s "{input.allreads_unfiltered}" ]; then
            echo "⚠️ Missing or empty input files for {wildcards.subset} — skipping"
            : > "{output}"
            exit 0
        fi

        usearch -otutab "{input.allreads_unfiltered}" \
                -zotus "{input.zotus}" \
                -otutabout "{output}" \
                -threads {threads} \
                -sample_delim "{params.sample_sep}"
        """


rule merge_abund_tables:
    input:
        lambda wc: expand(
            os.path.join(config["tmp_dir"], "05-abund_table", "abund_table_{subset}.tsv"),
            subset=get_dynamic_subsample_sizes(wc)
        )
    output:
        os.path.join(config["output_dir"], "07-abund_table", "abund_table_merged.tsv")
    log:
        os.path.join(config["log_dir"], "07-abund_table", "merge_abund_tables.log")
    message:
        "Merging abundance tables across all subsets"
    threads: 1
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p "$(dirname {output})"

        nfiles=$(echo {input} | wc -w)
        if [ "$nfiles" -eq 1 ]; then
            echo "⚠️ Only one abundance table found — copying instead of merging"
            cp {input} {output}
        else
            usearch -otutab_merge {input} -output "{output}"
        fi
        """



#rule rarefy_abund_table:
#    input:
#        os.path.join(config["output_dir"], "05-abund_table", "abund_table_merged.tsv")
#    output:
#        os.path.join(config["output_dir"], "05-abund_table", "abund_table_rarefied.tsv")
#    log:
#        os.path.join(config["log_dir"], "05-abund_table", "rarefy_abund_tables.log")
#    message:
#        "Rarefying merged abundance table"
#    threads: 1
#    container:
#        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
#    conda:
#        "../envs/snakemake_usearch.yml"
#    params:
#        rarefy_sample_size=config.get("rarefy_sample_size", 10000),
#    shell:
#        r"""
#        exec &> "{log}"
#        set -euxo pipefail
#
#        usearch -otutab_rare "{input}" -sample_size {params.rarefy_sample_size} -output "{output}"
#        """