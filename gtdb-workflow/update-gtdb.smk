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

# Dictionary for dynamic slurm batch allocations with correct resources
#PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 33], 4: ['med2', 33], 5: ['high2', 100]}
PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 33], 4: ['bmm', 33], 5: ['bmh', 100]}

wildcard_constraints:
    k = "\d{2}+",
    r = "\d[^-.]+",
    OR = "\w[^-.]+",

rule all:
    input:
        expand("{o}/gtdb-rs{r}-k{k}.zip", o=OUTDIR, r=RELEASES, k=KSIZES),
        expand("{o}/gtdb-rs{r}-k{k}.clean.zip", o=OUTDIR, r=RELEASES, k=KSIZES),
        expand("{o}/gtdb-rs{r}.lineages.csv", o=OUTDIR, r=RELEASES),
        expand("{o}/data/gtdb-{r}.clean-existing.csv", o=OUTDIR, r=RELEASES),
        expand("{o}/data/gtdb-{r}.clean-missing.csv", o=OUTDIR, r=RELEASES),

rule tax:
    input:
        expand("{o}/gtdb-rs{r}.lineages.csv", o=OUTDIR, r=RELEASES),

rule download_metadata:
    output:
        ar = '{o}/data/ar53_metadata_rs{r}.tsv',
        bac = '{o}/data/bac120_metadata_rs{r}.tsv',
    params:
        ar_path = lambda wildcards: f"https://data.gtdb.ecogenomic.org/releases/release{wildcards.r}/{wildcards.r}.{VERSION}",
        ar_file = lambda wildcards: f"ar53_metadata_r{wildcards.r}.tsv.gz",
        bac_path = lambda wildcards: f"https://data.gtdb.ecogenomic.org/releases/release{wildcards.r}/{wildcards.r}.{VERSION}",
        bac_file = lambda wildcards: f"bac120_metadata_r{wildcards.r}.tsv.gz",
    shell: """
        curl -L {params.bac_path}/{params.bac_file} > {wildcards.o}/data/{params.bac_file}
        curl -L {params.ar_path}/{params.ar_file} > {wildcards.o}/data/{params.ar_file}
        gzip -d {wildcards.o}/data/{params.bac_file} -c > {output.bac}
        gzip -d {wildcards.o}/data/{params.ar_file} -c > {output.ar}
    """

rule download_assembly_summary:
    output:
        good = '{o}/data/assembly_summary.bac.x.ar.txt',
        bad = '{o}/data/assembly_summary_historical.bac.x.ar.txt',
    resources:
        mem_mb = lambda wildcards, attempt: 2 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt > {wildcards.o}/data/assembly_summary.bacteria.txt
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt > {wildcards.o}/data/assembly_summary.archaea.txt

        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary_historical.txt > {wildcards.o}/data/assembly_summary_historical.bacteria.txt
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary_historical.txt > {wildcards.o}/data/assembly_summary_historical.archaea.txt

        cat {wildcards.o}/data/assembly_summary.bacteria.txt > {output.good}
        sed '1,2d' {wildcards.o}/data/assembly_summary.archaea.txt >> {output.good}
        cat {wildcards.o}/data/assembly_summary_historical.bacteria.txt > {output.bad}
        sed '1,2d' {wildcards.o}/data/assembly_summary_historical.archaea.txt >> {output.bad}

    """

rule get_ss_db:
    params:
        o_r = OLD_RELEASES
    output:
        dbs = temporary(f"gtdb-rs{OLD_RELEASES}-k{{k}}.zip"),
    resources:
        mem_mb = lambda wildcards, attempt: 2 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
            echo "Previous release: {params.o_r}"

            echo "Checking if gtdb-rs{params.o_r}-k{wildcards.k}.zip exists..."
            if [ -e /group/ctbrowngrp/sourmash-db/gtdb-rs{params.o_r}/gtdb-rs{params.o_r}-k{wildcards.k}.zip ]; then

                echo "gtdb-rs{params.o_r}-k{wildcards.k}.zip exists!"
                echo "Linking existing file to $(pwd)"

                ln -s /group/ctbrowngrp/sourmash-db/gtdb-rs{params.o_r}/gtdb-rs{params.o_r}-k{wildcards.k}.zip $(pwd)/{output.dbs}
            else

                echo "gtdb-rs{params.o_r}-k{wildcards.k}.zip does not exist!"
                echo "Downloading file to $(pwd)/{output.dbs}"

                curl -L https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs{params.o_r}/gtdb-rs{params.o_r}-k{wildcards.k}.zip > {output.dbs}
            fi
   """

rule make_taxonomy:
    input:
        script = 'scripts/make-gtdb-taxonomy.py',
        ar53_metadata='{o}/data/ar53_metadata_rs{r}.tsv',
        bac120_metadata='{o}/data/bac120_metadata_rs{r}.tsv',
    output:
        tax_csv = '{o}/gtdb-rs{r}.lineages.csv',
        reps_csv = '{o}/gtdb-rs{r}.lineages.reps.csv',
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    conda:
        "envs/sourmash.yaml",
    shell:
        """
        python {input.script} --metadata-files {input.ar53_metadata} {input.bac120_metadata} \
               -o {output.tax_csv} --reps-csv {output.reps_csv}
        """

rule collect_all:
    input:
        dbs = expand(f"gtdb-rs{OLD_RELEASES}-k{{k}}.zip", k = KSIZES)
    output:
        db = f"{{o}}/data/collect-mf.{OLD_RELEASES}.csv",
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=100,
        partition="bmh",
    shell: """
        sourmash sig collect {input.dbs} -F csv -o {output.db}
    """

rule picklist_check:
    input:
        dbs_manifest = f"{{o}}/data/collect-mf.{OLD_RELEASES}.csv",
        tax_picklist = '{o}/gtdb-rs{r}.lineages.csv',
    output:
        missing = "{o}/data/gtdb-{r}.missing.csv",
        manifest = "{o}/data/gtdb-{r}.existing.csv",
    params: 
        log = "logs/picklist_check/gtdb-{r}.picklist_check.log",
        first_k = lambda wildcards: KSIZES[0],
    conda: "envs/sourmash.yaml"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: 6 * 1024 * attempt,
        time= 10000,
        partition='high2',
    shell:
        """
        sourmash sig check -k {params.first_k} \
            --picklist {input.tax_picklist}:ident:ident \
            {input.dbs_manifest} --output-missing {output.missing} \
            --save-manifest {output.manifest} 2> {params.log}
        touch {output.missing}
        """

rule cleanse_existing:
    input:
        script = "scripts/update_sourmash_dbs.py",
        good = "{o}/data/assembly_summary.bac.x.ar.txt",
        bad = "{o}/data/assembly_summary_historical.bac.x.ar.txt",
        manifest = "{o}/data/gtdb-{r}.existing.csv",
    output:
        clean = "{o}/data/gtdb-{r}.clean-existing.csv",
        reversion = "{o}/data/gtdb-{r}.updated-versions-existing.csv",
        report = "{o}/data/report-existing.{r}.txt",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean} --updated-version {output.reversion} --report {output.report}
    """

rule cleanse_missing:
    input:
        script = "scripts/update_sourmash_dbs.py",
        good = "{o}/data/assembly_summary.bac.x.ar.txt",
        bad = "{o}/data/assembly_summary_historical.bac.x.ar.txt",
        manifest = "{o}/data/gtdb-{r}.missing.csv",
    output:
        clean = "{o}/data/gtdb-{r}.clean-missing.csv",
        all_links = "{o}/data/gtdb-{r}.all-missing-links.csv",
        report = "{o}/data/report-missing.{r}.txt",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean} --all-links {output.all_links} --report {output.report}
    """

rule gather_sketch_existing:
    input:
        missing = "{o}/data/gtdb-{r}.updated-versions-existing.csv",
    output:
        failed = "{o}/data/gtdb-{r}.updated-versions-existing.failures.csv",
        db = temporary("{o}/gtdb-{r}.updated-versions-existing.zip"),
    conda: "envs/directsketch.yaml"
    resources:
        mem_mb = 100 * 1024,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=100,
        partition="bmh",
    threads: 3
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        scale = config.get('scale_value'),
        log = "logs/gather_sketch_existing.{r}.log"
    shell:'''
        sourmash scripts gbsketch {input.missing} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled={params.scale},abund" -r 5 -g 2> {params.log}
    '''

rule gather_sketch_missing:
    input:
        missing = "{o}/data/gtdb-{r}.all-missing-links.csv",
    output:
        failed = "{o}/data/gtdb-{r}.all-missing-links.failures.csv",
        db = temporary("{o}/gtdb-{r}.all-missing-links.zip"),
    conda: "envs/directsketch.yaml"
    resources:
        mem_mb = 100 * 1024,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=100,
        partition="bmh",
    threads: 3
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        scale = config.get('scale_value'),
        log = "logs/gather_sketch_missing.{r}.log"
    shell:'''
        sourmash scripts gbsketch {input.missing} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled={params.scale},abund" -r 5 -g 2> {params.log}
    '''

rule extract_db:
    input:
        clean = "{o}/data/gtdb-{r}.clean-existing.csv",
        dbs = f"gtdb-rs{OLD_RELEASES}-k{{k}}.zip",
    output:
        dbs = "{o}/gtdb-rs{r}-k{k}.clean.zip",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = 80 * 1024,
        time = lambda wildcards, attempt: 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 24 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    shell:"""
        sourmash signature cat {input.dbs} --picklist {input.clean}:name:name -k {wildcards.k} -o {output.dbs}
    """

# Then need to combine everything into a single database
#        dbs = "{o}/gtdb-rs{r}-k{k}.clean.zip",
#        db = temporary("{o}/gtdb-{r}.clean-missing.zip"),
#        db = temporary("{o}/gtdb-{r}.updated-versions-missing.zip"),
#        db = temporary("{o}/gtdb-{r}.updated-versions-existing.zip"),

rule final_db:
    input:
        old_existing = "{o}/gtdb-rs{r}-k{k}.clean.zip",
        updated_existing = "{o}/gtdb-{r}.updated-versions-existing.zip",
        updated_missing = "{o}/gtdb-{r}.all-missing-links.zip",
    output:
        woohoo = "{o}/gtdb-rs{r}-k{k}.zip"
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = 80 * 1024,
        time = lambda wildcards, attempt: 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 24 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    shell:"""
        sourmash signature cat \
            {input.old_existing} {input.updated_existing} {input.updated_missing} \
            -k {wildcards.k} -o {output.woohoo}
    """
