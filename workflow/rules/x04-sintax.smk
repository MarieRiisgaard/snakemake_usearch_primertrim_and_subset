rule sintax:
    input:
        os.path.join(config["output_dir"], "zOTUs.fa"),
    output:
        os.path.join(config["output_dir"], "zOTUs.sintax"),
    log:
        os.path.join(config["log_dir"], "04-sintax", "sintax.log"),
    message:
        "Predicting taxonomy of zOTUs using SINTAX"
    resources:
        mem_mb=4096,
        runtime=60,
        cpus_per_task=config["max_threads"],
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    threads: config["max_threads"]
    params:
        db=config["db_sintax"],
    shell:
        """
        exec &> "{log}"
        set -euxo pipefail

        echo "Filtering zOTUs shorter than {params.minlen} bp before SINTAX"
        # --- Count before filtering ---
        total_before=$(grep -c "^>" "{input}" || true)

        # --- Clean and filter sequences ---
        awk -v minlen={params.minlen} '
            BEGIN { 
                header = ""; seq = "" 
            }
            /^>/ {
                # If there is a previous record, process it
                if (header != "") {
                    if (length(seq) >= minlen)
                        print header "\n" seq
                }
                # Start new record
                header = $0
                seq = ""
                next
            }
            {
                # Add sequence line (remove spaces, just in case)
                gsub(/[^A-Za-z]/, "", $0)
                seq = seq $0
            }
            END {
                # Process last record
                if (header != "" && length(seq) >= minlen)
                    print header "\n" seq
            }
        ' "{input}" > "{input}.filtered"

        # --- Run SINTAX classification ---
        usearch -sintax \
            "{input}" \
            -db "{params.db}" \
            -tabbedout "{output}" \
            -strand both \
            -sintax_cutoff 0.8 \
            -threads "{threads}"
        sort -V "{output}" -o "{output}"

        # error if output file is empty
        if [ ! -s "{output}" ]; then
            echo "output file {output} is empty, exiting!"
            exit 1
        fi

        echo "✅ SINTAX completed successfully for {wildcards.subset}"

        """
