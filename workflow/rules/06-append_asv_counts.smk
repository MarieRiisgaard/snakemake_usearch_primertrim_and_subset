############################################
# 06 – Append ASV counts and reference matches
############################################
rule append_asv_counts:
    input:
        summary=os.path.join(config["output_dir"], "subsample", "subsample_summary.tsv"),
        zotus=expand(
            os.path.join(config["output_dir"], "04-denoise", "{subset}", "zOTUs.fa.filtered"),
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
    
        mkdir -p $(dirname {output.summary_with_asvs})
    
        # --- If no summary file, stop gracefully ---
        if [ ! -s "{input.summary}" ]; then
            echo -e "Subset\tASVs\tASV_ref_matches\tAll_reads_ASVs\tAll_reads_ref_matches" > "{output.summary_with_asvs}"
            echo -e "Subset\tASV_ID\tReference_ID\tIdentity" > "{output.detailed_matches}"
            exit 0
        fi
    
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
            echo -e "${subset}\t${n_asv}" >> "$tmp_counts"
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
                    usearch -usearch_global "$zotus_file" \
                        -db "{params.reference_amplicons}" \
                        -id 0.95 \
                        -strand both \
                        -userout "${zotus_file}_refmatches.txt" \
                        -userfields query+target+id \
                        -quiet \
                        -maxaccepts 0 -maxrejects 0 \
                        -top_hits_only no || true
    
                    n_match=$(awk '$3==1.0' "${zotus_file}_refmatches.txt" | wc -l || echo 0)
    
                    if [ -s "${zotus_file}_refmatches.txt" ]; then
                        awk -v subset="$subset" '{print subset"\t"$1"\t"$2"\t"$3}' "${zotus_file}_refmatches.txt" >> "{output.detailed_matches}"
                    fi
                else
                    n_match=0
                fi
                echo -e "${subset}\t${n_match}" >> "$tmp_ref"
            done
        else
            echo "⚠️ No reference_amplicons specified or file missing — skipping identity checks."
            for zotus_file in {input.zotus}; do
                subset=$(basename "$(dirname "$zotus_file")")
                echo -e "${subset}\tNA" >> "$tmp_ref"
            done
        fi
    
        # --- Merge ASV counts + reference matches + all_reads data ---
        awk -v FS="\t" -v OFS="\t" -v cnt="$tmp_counts" -v ref="$tmp_ref" '
        BEGIN {
            while ((getline line < cnt) > 0) {
                if (line ~ /^Subset/) continue
                split(line, f, "\t")
                a[f[1]] = f[2]
            }
            close(cnt)
            while ((getline line < ref) > 0) {
                if (line ~ /^Subset/) continue
                split(line, f, "\t")
                r[f[1]] = f[2]
            }
            close(ref)
    
            all_asv = (("all_reads" in a) ? a["all_reads"] : 0)
            all_ref = (("all_reads" in r) ? r["all_reads"] : "NA")
        }
        NR == 1 {
            print $0, "ASVs", "ASV_ref_matches", "All_reads_ASVs", "All_reads_ref_matches"
            next
        }
        /^barcode/ {
            subset = "sample_size_" $2
            asv = (subset in a ? a[subset] : 0)
            refv = (subset in r ? r[subset] : (("all_reads" in r) ? r["all_reads"] : "NA"))
            print $0, asv, refv, all_asv, all_ref
            next
        }
        {print}
        ' "{input.summary}" > "{output.summary_with_asvs}"
    
        rm "$tmp_counts" "$tmp_ref"
        echo "✅ append_asv_counts completed successfully."
        """
    