rule append_asv_counts:
    input:
        summary=os.path.join(config["output_dir"], "subsample", "subsample_summary.tsv"),
        zotus=expand(
            os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa"),
            subset=subset_dirs
        )
    output:
        summary_with_asvs=os.path.join(config["output_dir"], "subsample_summary_with_ASVs.tsv"),
        detailed_matches=os.path.join(config["output_dir"], "subsample_summary_ASV_matches.tsv")
    log:
        os.path.join(config["log_dir"], "06-append_asv_counts.log")
    params:
        reference_amplicons=config.get("reference_amplicons", "")
    message:
        "Appending ASV counts and generating reference match report"
    container:
        "docker://ghcr.io/kasperskytte/snakemake_usearch:main"
    conda:
        "../envs/snakemake_usearch.yml"
    shell:
        r"""
        exec &> "{log}"
        set -euxo pipefail

        # --- Collect ASV counts per subset ---
        tmp_counts=$(mktemp)
        echo -e "Subset\tASVs" > "$tmp_counts"
        for zotus_file in {input.zotus}; do
            subset=$(basename "$(dirname "$zotus_file")")
            if [ -s "$zotus_file" ]; then
                n_asv=$(grep -c '^>' "$zotus_file" || echo 0)
            else
                n_asv=0
            fi
            echo -e "${{subset}}\t${{n_asv}}" >> "$tmp_counts"
        done

        # --- Prepare reference match report ---
        tmp_ref=$(mktemp)
        echo -e "Subset\tASV_ref_matches" > "$tmp_ref"
        echo -e "Subset\tASV_ID\tReference_ID\tIdentity" > "{output.detailed_matches}"

        if [ -n "{params.reference_amplicons}" ] && [ -s "{params.reference_amplicons}" ]; then
            echo "Reference file found: {params.reference_amplicons}"
            for zotus_file in {input.zotus}; do
                subset=$(basename "$(dirname "$zotus_file")")
                if [ -s "$zotus_file" ]; then
                    # Usearch full search for all hits ≥0.9 identity
                    usearch -usearch_global "$zotus_file" \
                        -db "{params.reference_amplicons}" \
                        -id 0.95 \
                        -strand both \
                        -userout "${{zotus_file}}_refmatches.txt" \
                        -userfields query+target+id \
                        -quiet \
                        -maxaccepts 0 -maxrejects 0 \
                        -top_hits_only no

                    # Count full 100% matches
                    n_match=$(awk '$3==1.0' "${{zotus_file}}_refmatches.txt" | wc -l || echo 0)

                    # Append all hits (full + partial) to the detailed output
                    if [ -s "${{zotus_file}}_refmatches.txt" ]; then
                        awk -v subset="$subset" '{{print subset"\t"$1"\t"$2"\t"$3}}' "${{zotus_file}}_refmatches.txt" >> "{output.detailed_matches}"
                    fi
                else
                    n_match=0
                fi
                echo -e "${{subset}}\t${{n_match}}" >> "$tmp_ref"
            done
        else
            echo "⚠️  No reference_amplicons specified or file missing — skipping identity checks."
            for zotus_file in {input.zotus}; do
                subset=$(basename "$(dirname "$zotus_file")")
                echo -e "${{subset}}\tNA" >> "$tmp_ref"
            done
        fi

        # --- Merge ASV counts + reference matches into main summary ---
        awk '
            BEGIN {{FS=OFS="\t"}}
            FNR==NR {{asv[$1]=$2; next}}
            NR==FNR+1 {{ref[$1]=$2; next}}
            NR==1 {{print $0, "ASVs", "ASV_ref_matches"; next}}
            /^barcode/ {{
                subset="sample_size_"$2
                print $0, (asv[subset]?asv[subset]:0), (ref[subset]?ref[subset]:"NA")
                next
            }}
        ' "$tmp_counts" "$tmp_ref" "{input.summary}" > "{output.summary_with_asvs}"

        # Add all_reads ASV/ref count for subsets missing a size
        all_asv=$(awk 'NR>1 && $1=="all_reads" {{print $2}}' "$tmp_counts" || echo 0)
        all_ref=$(awk 'NR>1 && $1=="all_reads" {{print $2}}' "$tmp_ref" || echo "NA")
        awk -v asv="$all_asv" -v ref="$all_ref" 'BEGIN{{FS=OFS="\t"}} NR==1{{print $0; next}} !/\tsample_size_/{{print $0,asv,ref; next}} {{print}}' "{output.summary_with_asvs}" > "{output.summary_with_asvs}.tmp" && mv "{output.summary_with_asvs}.tmp" "{output.summary_with_asvs}"

        rm "$tmp_counts" "$tmp_ref"
        """
