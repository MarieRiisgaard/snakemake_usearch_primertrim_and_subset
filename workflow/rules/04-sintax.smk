rule sintax:
    input:
        os.path.join(config["output_dir"], "zOTUs.fa")
    output:
        os.path.join(config["output_dir"], "zOTUs.sintax")
    log:
        os.path.join(config["log_dir"], "04-sintax", "sintax.log")
    message:
        "Predicting taxonomy of zOTUs using SINTAX",
    resources:
        mem_mb=4096,
        runtime=60,
        cpus_per_task=config["max_threads"]
    threads: config["max_threads"]
    conda:
        "../envs/env.yml"
    params:
      db=config["db_sintax"]
    shell:
      """
        exec &> "{log}"
        set -euxo pipefail
            
        usearch11 -sintax \
          "{input}" \
          -db "{params.db}" \
          -tabbedout "{output}" \
          -strand both \
          -sintax_cutoff 0.8 \
          -threads "{threads}"
        sort -V "{output}" -o "{output}"
      """
