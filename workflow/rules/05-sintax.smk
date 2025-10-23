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
        db=config["db_sintax"],
        minlen=config.get("minlen", 100)
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail
    
        echo "Filtering zOTUs shorter than {params.minlen} bp before SINTAX"
    
        total_before=$(grep -c "^>" "{input}" || true)
    
        awk -v minlen={params.minlen} '
            BEGIN {{
                header = ""; seq = ""
            }}
            /^>/ {{
                if (header != "") {{
                    if (length(seq) >= minlen)
                        print header "\\n" seq
                }}
                header = $0
                seq = ""
                next
            }}
            {{
                gsub(/[^A-Za-z]/, "", $0)
                seq = seq $0
            }}
            END {{
                if (header != "" && length(seq) >= minlen)
                    print header "\\n" seq
            }}
        ' "{input}" > "{input}.filtered"
    
        usearch -sintax \
            "{input}.filtered" \
            -db "{params.db}" \
            -tabbedout "{output}" \
            -strand both \
            -sintax_cutoff 0.8 \
            -threads "{threads}"
    
        sort -V "{output}" -o "{output}"
    
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi
    
        echo "✅ SINTAX completed successfully for {wildcards.subset}"
        """