# Configuration
The configuration file `config.yaml` is used to set various options used throughout the workflow.

| Option | Default value | Description |
| ---: | :---: | :---: |
| input_dir | `"data/"` | The input folder is expected to contain a subfolder for each sampleID/barcode, in which all fastq files will be concatenated, and the subfolder names used as sample IDs downstream. For nanopore this is usually the "fastq_pass" folder with demultiplexed reads. |
| output_dir | `"output"` | Folder for the results. |
| tmp_dir | `"tmp"` | Folder for temporary files, which are deleted by default after a succesful run. |
| log_dir | `"logs"` | Folder for logs for each rule. |
| db_sintax | `"/databases/midas/MiDAS5.3_20240320/FLASVs_w_sintax.fa"` | Path to the taxonomic reference database used to classify the ASVs/zOTUs in SINTAX format. |
| filtlong_args | `"--min_length 200 --min_mean_q 90"` | Arguments for the filtlong command used for pre-filtering. To skip filtering altogether set to `"--min_length 1"`. |
| max_threads | `32` | Max number of threads to use for any individual rule. |
| sample_sep | `"_"` | Separator used for the `usearch -otutab` and `fastx_relabel` commands. |
| primers | `AGRGTTYGATYMTGGCTCAG...GTTTGGCACCTCGATGTCG"` | Primer pair used. Passed on as-is to the `cutadapt` command. This is required for trimming and orienting reads correctly. |
| derep_minsize | `2` | Minimum abundance of each read. Only impacts ASV/zOTU generation, not abundance estimation, which will be performed against the raw unfiltered reads. |
| unoise_minsize | `8` | Increase this proportionally with platform error-rate to avoid false-positive de-novo ASVs/zOTUs. Never set to lower than 2, singletons cannot be trusted. |
| rarefy | `2000` | Rarefy abundance table to an equal sample size. Both a rarefied and an unrarefied abundance table will be generated. |

Have a look in the `.test` directory for minimal example files.
