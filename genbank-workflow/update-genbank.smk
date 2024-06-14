###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s update-genbank.smk -j 15 --use-conda --rerun-incomplete --resources allowed_jobs=100
#
# On HPC, use 10 cpus and ~50gb. Request 2 days?
###

configfile: "config/update-genbank.yaml"

DATE = config.get("update_to_date")

if not DATE:
    import time
    DATE = [time.strftime("%Y%m%d")]

#f"/group/ctbrowngrp/sourmash-db/genbank-{DATE}"
outdir = [
    f"{config.get('output_directory') if config.get('output_directory') is not None else '..'}/genbank-{date}"
    for date in DATE]

DOMAINS = config.get('domains')

KSIZES = config.get('k_values')

email = config.get('email')

if email:
    onsuccess:
        print("\nWorkflow finished without error\n")
        shell("mail -s 'Workflow finished without error' {email} < {log}")

    onerror:
        print("\nAn error occurred\n")
        shell("mail -s 'an error occurred' {email} < {log}")

OLD_DATES, = config["update_from_date"]

wildcard_constraints:
    k = "\d{2}+",
    D = "\w[^-.]+",
    d = "\d+",

rule all:
    input:
        expand("{o}/genbank-{d}-{D}-k{k}.zip", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),
        expand("{o}/lineages.{D}.csv", o=outdir, D=DOMAINS),

rule build_genbank:
    input:
        expand("{o}/genbank-{d}-{D}-k{k}.zip", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),

rule clean_genbank:
    input:
        expand("{o}/genbank-{d}-{D}-k{k}.clean.zip", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),

rule missing_genbank:
    input:
        expand("{o}/genbank-{d}-{D}-k{k}.missing.zip", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),

rule check_genbank:
    input:
        expand("{o}/data/genbank-{d}-{D}-k{k}.zip.check", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),

rule tax_genbank:
    input:
        expand("{o}/lineages.{D}.csv", o=outdir, D=DOMAINS),

rule download_assembly_summary:
    output:
        good = '{o}/data/assembly_summary.{D}.txt',
        bad = '{o}/data/assembly_summary_historical.{D}.txt',
    benchmark:
        'benchmarks/download_assembly_summary.{o}_{D}.tsv'
    shell: """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary.txt > {output.good}
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary_historical.txt > {output.bad}
    """

rule get_ss_db:
    params:
        old_date = OLD_DATES
    output:
        dbs = temporary(f"genbank-{OLD_DATES}-{{D}}-k{{k}}.zip"),
    conda: "envs/sourmash.yaml"
    shell: """
            echo "Old date: {params.old_date}"

            echo "Checking if genbank-{params.old_date}-{wildcards.D}-k{wildcards.k}.zip exists..."
            if [ -e /group/ctbrowngrp/sourmash-db/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.D}-k{wildcards.k}.zip ]; then

                echo "genbank-{params.old_date}-{wildcards.D}-k{wildcards.k}.zip exists!"
                echo "Linking existing file to $(pwd)"

                ln -s /group/ctbrowngrp/sourmash-db/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.D}-k{wildcards.k}.zip $(pwd)/{output.dbs}
            else

                echo "genbank-{params.old_date}-{wildcards.D}-k{wildcards.k}.zip does not exist!"
                echo "Downloading file to $(pwd)/{output.dbs}"

                curl -L https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.D}-k{wildcards.k}.zip > {output.dbs}
            fi
   """

# Does this need to be updated because some scripts are not on main branches?
rule download_update_sourmash_dbs:
    output: "scripts/update_sourmash_dbs.py"
    shell: """
        curl -L "https://raw.githubusercontent.com/ccbaumler/2024-database-create/manifests/workflow/scripts/update_sourmash_dbs.py" > {output}
    """

rule manifest_manifest:
    input:
        dbs = f"genbank-{OLD_DATES}-{{D}}-k{{k}}.zip",
    output: "{o}/data/mf.{d}-{D}-k{k}.csv",
    conda: "envs/sourmash.yaml",
    params:
        old_date=OLD_DATES,
    benchmark:
        'benchmarks/mf.{o}-{d}-{D}-k{k}.tsv'
    shell: """
        for db_file in {input}; do
            if [[ $db_file == *{wildcards.D}*{wildcards.k}* ]]; then
                sourmash sig manifest -o {output} $db_file --no-rebuild
            fi
        done
    """

rule cleanse_manifest:
    input:
        script = "scripts/update_sourmash_dbs.py",
        good = "{o}/data/assembly_summary.{D}.txt",
        bad = "{o}/data/assembly_summary_historical.{D}.txt",
        manifest = "{o}/data/mf.{d}-{D}-k{k}.csv",
    output:
        clean = "{o}/data/mf-clean.{d}-{D}-k{k}.csv",
        reversion = "{o}/data/updated-versions.{d}-{D}-k{k}.csv",
        report = "{o}/data/report.{d}-{D}-k{k}.txt",
        missing = "{o}/data/missing-genomes.{d}-{D}-k{k}.csv",
    conda: "envs/sourmash.yaml",
    benchmark:
        'benchmarks/cleanse_manifest.{o}-{d}-{k}-{D}.tsv'
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean} --updated-version {output.reversion} --report {output.report} --missing-genomes {output.missing}
    """

rule picklist_clean_db:
    input:
        clean = "{o}/data/mf-clean.{d}-{D}-k{k}.csv",
        dbs = f"genbank-{OLD_DATES}-{{D}}-k{{k}}.zip",
    output:
        woohoo = temporary("{o}/genbank-{d}-{D}-k{k}.clean.zip"),
    conda: "envs/sourmash.yaml",
    params:
        old_date = OLD_DATES
    benchmark:
        'benchmarks/picklist_clean_db.{o}-{d}-{D}-k{k}.tsv'
    shell:'''
        if [ ! -e {output.woohoo} ]; then
            echo "Cleaning {input.dbs}..."
            sourmash sig extract --picklist {input.clean}::manifest {input.dbs} -o {output.woohoo}
            echo "{input.dbs} cleaned and stored as {output.woohoo}"
        fi
    '''

rule check_txt_reversioned:
    input:
        reversion = expand("{o}/data/updated-versions.{d}-{D}-k{k}.csv", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),
        script = "scripts/check_txt_files.py",
    output:
        solo = "{o}/data/update.{d}-{D}.csv",
    benchmark:
        'benchmarks/check_txt_reversioned.{o}-{d}-{D}.tsv'
    shell: """
        files=""
        for file in "{wildcards.o}/data/updated-versions.*{wildcards.d}-{wildcards.D}-k*.csv"; do
            files+=" $file"
        done
        {input.script} $files -o {output.solo}
    """

rule gather_sketch_reversioned:
    input:
        reversion = "{o}/data/update.{d}-{D}.csv",
    output:
        failed = "{o}/data/update.{d}-{D}.failures.csv",
        db = temporary("{o}/genbank-{d}-{D}.rever.zip"),
    conda: "envs/directsketch.yaml"
    benchmark:
        'benchmarks/gather_sketch_reversioned.{o}-{d}-{D}.tsv'
    threads: 3
    resources:
        allowed_jobs = 100
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
    log:
        "logs/gather_sketch_reversioned.{o}_{d}_{D}.log"
    shell:'''
        sourmash scripts gbsketch {input.reversion} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled=1000,abund" -r 5 -g 2> {log}
    '''

rule cat_to_clean_reversioned:
    input:
        dir = "{o}/genbank-{d}-{D}.rever.zip",
        db = "{o}/genbank-{d}-{D}-k{k}.clean.zip",
        missing = "{o}/data/update.{d}-{D}.failures.csv",
    output:
        woohoo = temporary("{o}/genbank-{d}-{D}-k{k}.update.zip"),
    conda: "envs/sourmash.yaml"
    benchmark:
        'benchmarks/cat_to_clean_reversioned.{o}-{d}-{D}-k{k}.tsv'
    shell: """
        sourmash sig cat {input.dir} {input.db} -k {wildcards.k} -o {output.woohoo}
    """

rule check_txt_missing:
    input:
        missing = expand("{o}/data/missing-genomes.{d}-{D}-k{k}.csv", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),
        script = "scripts/check_txt_files.py",
    output:
        solo = "{o}/data/missing.{d}-{D}.csv",
    benchmark:
        'benchmarks/check_txt_missing.{o}-{d}-{D}.tsv'
    shell: """
        files=""
        for file in "{wildcards.o}/data/missing-genomes.*{wildcards.d}-{wildcards.D}-k*.csv"; do
            files+=" $file"
        done
        {input.script} $files -o {output.solo}
    """

rule gather_sketch_missing:
    input:
         missing = "{o}/data/missing.{d}-{D}.csv",
    output:
        failed = "{o}/data/missing.{d}-{D}.failures.csv",
        db = temporary("{o}/genbank-{d}-{D}.miss.zip"),
    conda: "envs/directsketch.yaml"
    benchmark:
        'benchmarks/gather_sketch_missing.{o}-{d}-{D}.tsv'
    threads: 3
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
    log:
        "logs/gather_sketch_missing.{o}_{d}_{D}.log"
    resources:
        allowed_jobs = 100
    shell:'''
        sourmash scripts gbsketch {input.missing} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled=1000,abund" -r 5 -g 2> {log}
    '''

checkpoint cat_to_clean_missing:
    input:
        dir = "{o}/genbank-{d}-{D}.miss.zip",
        db = "{o}/genbank-{d}-{D}-k{k}.clean.zip",
        missing = "{o}/data/missing.{d}-{D}.failures.csv",
    output:
        woohoo = protected("{o}/genbank-{d}-{D}-k{k}.zip"),
    benchmark:
        'benchmarks/cat_to_clean_missing.{o}-{d}-{D}-k{k}.tsv',
    conda: "envs/sourmash.yaml"
    resources:
        allowed_jobs = 100
    shell: """
        sourmash sig cat {input.dir} {input.db} -k {wildcards.k} -o {output.woohoo}
    """

### create a report with sourmash sig summarize for the databases... and sourmash compare(?)

# taxonomy rules, from https://github.com/ctb/2022-assembly-summary-to-lineages
rule download_ncbi_utils:
    output: "scripts/ncbi_taxdump_utils.py"
    shell:
        "curl -L https://raw.githubusercontent.com/ctb/2022-assembly-summary-to-lineages/main/ncbi_taxdump_utils.py > {output}"

rule download_taxscript:
    output: "scripts/make-lineage-csv.py"
    shell:
        "curl -L https://raw.githubusercontent.com/bluegenes/2022-assembly-summary-to-lineages/virus-tax/make-lineage-csv.py > {output}"

rule download_taxdump: # may need to restart this a couple times
    output:
        "taxdump/nodes.dmp",
        "taxdump/names.dmp"
    shell:
        "curl -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | (mkdir -p taxdump && cd taxdump && tar xzvf -)"

rule make_lineage_csv:
    input:
        "{o}/data/assembly_summary.{D}.txt",
        "taxdump/nodes.dmp",
        "taxdump/names.dmp",
        "scripts/make-lineage-csv.py",
        "scripts/ncbi_taxdump_utils.py",
    output:
        "{o}/lineages.{D}.csv"
    params:
        ictv_cmd = lambda w: " --ictv " if 'viral' in w.D else '',
    benchmark:
        'benchmarks/make_lineage_csv.{o}-{D}.tsv'
    shell:
        "python scripts/make-lineage-csv.py taxdump/{{nodes.dmp,names.dmp}} {input[0]} -o {output} {params.ictv_cmd}"
