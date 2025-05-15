
  scriptMessage "Filtering ASVs that are <60% similar to reference reads..."
  if [ -s "$prefilterdb" ]
  then
    usearch11 -usearch_global \
      "${tempdir}/zOTUs.R1.fa" \
      -db "$prefilterdb" \
      -strand both \
      -id 0.6 \
      -maxaccepts 1 \
      -maxrejects 8 \
      -matched "${tempdir}/prefilt_out.fa" \
      -threads "$maxthreads" \
      -quiet
    mv "${tempdir}/prefilt_out.fa" "${tempdir}/zOTUs.R1.fa"
  else
    echo "Could not find prefilter reference database, continuing without prefiltering..."
  fi

  scriptMessage "Searching ASVs against already known ASVs (exact match) and renaming accordingly..."
  if [ -s "$asvdb" ]
  then
    usearch11 -search_exact \
      "${tempdir}/zOTUs.R1.fa" \
      -db "$asvdb" \
      -maxaccepts 0 \
      -maxrejects 0 \
      -strand both \
      -dbmatched "${output}/ASVs.R1.fa" \
      -notmatched "${tempdir}/ASVs_nohits.R1.fa" \
      -threads "$maxthreads" \
      -quiet
    usearch11 -fastx_relabel \
      "${tempdir}/ASVs_nohits.R1.fa" \
      -prefix "${newASVprefix}" \
      -fastaout "${tempdir}/ASVs_nohits_renamed.R1.fa" \
      -quiet
    #combine hits with nohits
    cat "${tempdir}/ASVs_nohits_renamed.R1.fa" >> "${output}/ASVs.R1.fa"
  else
    echo "Could not find ASV database, continuing without renaming ASVs..."
    sed 's/Zotu/ASV/g' "${tempdir}/zOTUs.R1.fa" > "${output}/ASVs.R1.fa"
  fi
