###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s update-genbank.smk -j 15 --use-conda --rerun-incomplete --resources allowed_jobs=100
#
# On HPC, use 10 cpus and ~50gb. Request 2 days?
###

configfile: "config/update-gtdb.yaml"

RELEASES = config.get("update_to_release")
OLD_RELEASES, = config["update_from_release"]
VERSION = 0 if config.get("release_version") is None else config.get("release_version")

OUTDIR = [
    f"{config.get('output_directory') if config.get('output_directory') is not None else '..'}/gtdb-{release}"
    for release in RELEASES]

KSIZES = config.get('k_values')
EMAIL = config.get('email')

if EMAIL:
    onsuccess:
        print("\nWorkflow finished without error\n")
        shell("mail -s 'Workflow finished without error' {EMAIL} < {log}")

    onerror:
        print("\nAn error occurred\n")
        shell("mail -s 'an error occurred' {EMAIL} < {log}")

wildcard_constraints:
    k = "\d{2}+",
    R = "\w[^-.]+",
    OR = "\w[^-.]+",

rule all:
    input:
#        expand("{o}/gtdb-rs{r}-k{k}.zip", o=OUTDIR, r=RELEASES, k=KSIZES),
        expand("{o}/gtdb-rs{r}.lineages.csv", o=OUTDIR, r=RELEASES),

rule tax:
    input:
        expand("{o}/gtdb-rs{r}.lineages.csv", o=OUTDIR, r=RELEASES),

rule download_metadata:
    output:
        ar = '{o}/data/ar53_metadata_rs{r}.tsv',
        bac = '{o}/data/bac120_metadata_rs{r}.tsv',
    params:
        ar = lambda wildcards: f"https://data.gtdb.ecogenomic.org/releases/release{wildcards.r}/{wildcards.r}.{VERSION}/ar53_metadata_r{wildcards.r}.tsv",
        bac = lambda wildcards: f"https://data.gtdb.ecogenomic.org/releases/release{wildcards.r}/{wildcards.r}.{VERSION}/bac120_metadata_r{wildcards.r}.tsv",
    shell: """
        curl -L {params.bac} > {output.bac}
        curl -L {params.ar} > {output.ar}
    """


rule make_taxonomy:
    input:
        ar53_metadata='{o}/data/ar53_metadata_rs{r}.tsv',
        bac120_metadata='{o}/data/bac120_metadata_rs{r}.tsv',
    output:
        tax_csv = '{o}/gtdb-rs{r}.lineages.csv',
        reps_csv = '{o}/gtdb-rs{r}.lineages.reps.csv',
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 240,
        partition='high2',
    shell:
        """
        python make-gtdb-taxonomy.py --metadata-files {input.ar53_metadata} {input.bac120_metadata} \
               -o {output.tax_csv} --reps-csv {output.reps_csv}
        """
