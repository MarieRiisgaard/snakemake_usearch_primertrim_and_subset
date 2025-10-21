rule subsample_reads:
    input:
        os.path.join(
            config["tmp_dir"],
            "02-denoise",
            "all_samples_filtered_renamed_oriented_trimmed.fastq"
        )
    output:
        summary = os.path.join(config["output_dir"], "subsample", "{sample}", "subsample_summary.tsv")
    params:
        sizes = config.get("subsample_sizes", [1000, 2000, 5000, 10000])
    log:
        os.path.join(config["log_dir"], "03-subsample", "{sample}.log")
    conda:
        "../envs/snakemake_usearch.yml"
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    threads: 4
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p {config[output_dir]}/subsample/{wildcards.sample}

        # Concatenate all FASTQs
        concat_file={config[output_dir]}/subsample/{wildcards.sample}/{wildcards.sample}_all.fastq
        gunzip -cdfq {input} > $concat_file

        nreads=$(grep -c '^+$' $concat_file || true)
        echo -e "Sample\tSubsample_size\tReads_total\tReads_subsampled" > {output.summary}

        for size in {params.sizes}; do
            if [ "$nreads" -lt "$size" ]; then
                echo -e "{wildcards.sample}\t$size\t$nreads\t0" >> {output.summary}
                continue
            fi

            subsample_dir={config[output_dir]}/subsample/sample_size_${{size}}/{wildcards.sample}
            mkdir -p $subsample_dir
            subsample_file=$subsample_dir/{wildcards.sample}_subsampled_${{size}}.fastq

            usearch11 -fastx_subsample $concat_file \
                -sample_size $size \
                -fastqout $subsample_file \
                -randseed 42 \
                -quiet

            reads_after=$(grep -c '^+$' $subsample_file || true)
            echo -e "{wildcards.sample}\t$size\t$nreads\t$reads_after" >> {output.summary}
        done
        """
