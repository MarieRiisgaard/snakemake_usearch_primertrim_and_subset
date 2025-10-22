###############################################################
# 03-subsample.smk — unified structure with one summary file
###############################################################

# Subsample each barcode individually for each size and also copy full reads
rule subsample_reads:
    input:
        os.path.join(config["tmp_dir"], "02-denoise", "{sample}", "{sample}_trimmed.fastq")
    output:
        subsampled=os.path.join(
            config["output_dir"],
            "subsample",
            "sample_size_{size}",
            "{sample}_subsampled_{size}.fastq"
        ),
        all_reads=os.path.join(
            config["output_dir"],
            "subsample",
            "all_reads",
            "{sample}_trimmed_{size}.fastq"   # <-- add size for Snakemake consistency
        )
    log:
        os.path.join(config["log_dir"], "03-subsample", "{sample}_{size}.log")
    params:
        sizes=config.get("subsample_sizes", [10000, 20000, 50000, 100000, 200000, 400000])
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: 2
    message:
        "Subsampling {wildcards.sample} to size {wildcards.size}"
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        mkdir -p {config[output_dir]}/subsample/all_reads
        cp {input} {output.all_reads}

        input_file={input}
        nreads=$(grep -c '^+$' "$input_file" || true)

        size={wildcards.size}
        if [ "$nreads" -lt "$size" ]; then
            echo "⚠️  {wildcards.sample}: not enough reads ($nreads) for size $size, skipping..."
            exit 0
        fi

        subsample_dir={config[output_dir]}/subsample/sample_size_${{size}}
        mkdir -p $subsample_dir

        usearch -fastx_subsample "$input_file" \
            -sample_size $size \
            -fastqout {output.subsampled} \
            -randseed 42 \
            -quiet
        """


###############################################################
# Merge all subsampling results into one summary file
###############################################################

rule merge_subsample_summaries:
    input:
        expand(
            os.path.join(config["output_dir"], "subsample", "sample_size_{size}", "{sample}_subsampled_{size}.fastq"),
            sample=sample_dirs,
            size=config.get("subsample_sizes", [1000, 2000, 5000, 10000])
        )
    output:
        os.path.join(config["output_dir"], "subsample", "subsample_summary.tsv")
    params:
        samples=sample_dirs,
        sizes=config.get("subsample_sizes", [1000, 2000, 5000, 10000]),
        sample_list=" ".join(sample_dirs)
    log:
        os.path.join(config["log_dir"], "03-subsample", "merge_summary.log")
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        echo -e "Sample\tSubsample_size\tReads_total\tReads_trimmed\tReads_subsampled" > {output}
        echo "Samples: {params.sample_list}"
        echo "Sizes: {params.sizes}"

        for sample in {params.sample_list}; do
            # ----- Total reads BEFORE primer trimming -----
            reads_total_file={config[tmp_dir]}/totalreads/${{sample}}_totalreads.csv
            if [ -s "$reads_total_file" ]; then
                reads_total=$(awk -F',' 'NR==2 {{print $2}}' "$reads_total_file" || echo 0)
            else
                echo "⚠️  No total reads file for $sample"
                reads_total=0
            fi

            # ----- Reads AFTER primer trimming -----
            reads_trimmed_file={config[tmp_dir]}/totalreads_trimmed/${{sample}}_totaltrimmedreads.csv
            if [ -s "$reads_trimmed_file" ]; then
                reads_trimmed=$(awk -F',' 'NR==2 {{print $2}}' "$reads_trimmed_file" || echo 0)
            else
                echo "⚠️  No trimmed reads file for $sample"
                reads_trimmed=0
            fi

            # ----- Subsampled reads for each size -----
            for size in {params.sizes}; do
                subsample_file={config[output_dir]}/subsample/sample_size_${{size}}/${{sample}}_subsampled_${{size}}.fastq
                if [ -s "$subsample_file" ]; then
                    reads_subsampled=$(grep -c '^+$' "$subsample_file" || true)
                else
                    reads_subsampled=0
                fi
                echo -e "${{sample}}\t${{size}}\t${{reads_total}}\t${{reads_trimmed}}\t${{reads_subsampled}}" >> {output}
            done
        done
        """

###############################################################



#rule subsample_reads:
#    input:
#        os.path.join(
#            config["tmp_dir"],
#            "02-denoise",
#            "all_samples_filtered_renamed_oriented_trimmed.fastq"
#        )
#    output:
#        summary = os.path.join(config["output_dir"], "subsample", "{sample}", "subsample_summary.tsv")
#    params:
#        sizes = config.get("subsample_sizes", [1000, 2000, 5000, 10000])
#    log:
#        os.path.join(config["log_dir"], "03-subsample", "{sample}.log")
#    conda:
#        "../envs/snakemake_usearch.yml"
#    container:
#        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
#    threads: 4
#    shell:
#        r"""
#        exec &> "{log}"
#        set -euxo pipefail
#
#        mkdir -p {config[output_dir]}/subsample/{wildcards.sample}
#
#        # Concatenate all FASTQs
#        concat_file={config[output_dir]}/subsample/{wildcards.sample}/{wildcards.sample}_all.fastq
#        gunzip -cdfq {input} > $concat_file
#
#        nreads=$(grep -c '^+$' $concat_file || true)
#        echo -e "Sample\tSubsample_size\tReads_total\tReads_subsampled" > {output.summary}
#
#        for size in {params.sizes}; do
#            if [ "$nreads" -lt "$size" ]; then
#                echo -e "{wildcards.sample}\t$size\t$nreads\t0" >> {output.summary}
#                continue
#            fi
#
#            subsample_dir={config[output_dir]}/subsample/sample_size_${{size}}/{wildcards.sample}
#            mkdir -p $subsample_dir
#            subsample_file=$subsample_dir/{wildcards.sample}_subsampled_${{size}}.fastq
#
#            usearch -fastx_subsample $concat_file \
#                -sample_size $size \
#                -fastqout $subsample_file \
#                -randseed 42 \
#                -quiet
#
#            reads_after=$(grep -c '^+$' $subsample_file || true)
#            echo -e "{wildcards.sample}\t$size\t$nreads\t$reads_after" >> {output.summary}
#        done
#        """
#