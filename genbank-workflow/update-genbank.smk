###
# This workflow will create a genbanksize database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s update-genbank.smksize -j 15 --use-conda --rerun-incomplete --resources allowed_jobs=100
#
# On HPC, use 10 cpus and ~50gb. Request 2 days?
###

import os

NCBI_API_KEY = os.environ.get("NCBI_API_KEY")

configfile: "config/update-genbank.yaml"

DATE = config.get("update_to_date")

if not DATE:
    import time
    DATE = [time.strftime("%Y%m%d")]

INDIR = config.get("input_directory")
outdir = [
    f"{config.get('output_directory') if config.get('output_directory') is not None else '..'}/genbank-{date}"
    for date in DATE]

DOMAINS = config.get('domains')

KSIZES = config.get('k_values')
print(KSIZES)
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
    ksize = "\\d{2}+",
    dom = "\w[^-.]+",
    d = "\\d+",

# Dictionary for dynamic slurm batch allocations with correct resources
#PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 33], 4: ['med2', 33], 5: ['high2', 100]}
PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 33], 4: ['bmm', 33], 5: ['bmh', 100]}

# Create a file dictionary for normal and test runs for rule cheat_mainfest
def getInputFilesForManifest(wildcards):
    files = dict()
    if config.get('output_directory') == 'test':
        files["good"] = f"{wildcards.o}/data/sub_assembly_summary.{wildcards.dom}.txt"
        files["bad"] = f"{wildcards.o}/data/sub_assembly_summary_historical.{wildcards.dom}.txt"
    else:
        files["good"] = f"{wildcards.o}/data/assembly_summary.{wildcards.dom}.txt"
        files["bad"] = f"{wildcards.o}/data/assembly_summary_historical.{wildcards.dom}.txt"
    return files

# Create a file dictionary for normal and test runs for rule cheat_mainfest
def createOldSingleManifest(wildcards):
    file = dict()

    first_ksize = KSIZES[0] if isinstance(KSIZES, (list, tuple)) else KSIZES

    file["original_dbs"] = f"genbank-{OLD_DATES}-{wildcards.dom}-k{first_ksize}.zip"

    return file

def createNewSingleManifest(wildcards):
    file = dict()

    first_ksize = KSIZES[0] if isinstance(KSIZES, (list, tuple)) else KSIZES

    file["new_dbs"] = f"{wildcards.o}/genbank-{wildcards.d}-{wildcards.dom}-k{first_ksize}.zip"

    return file


#### Psuedo Rules ####


rule all:
    input:
        expand("{o}/genbank-{d}-{dom}-k{ksize}.zip", o=outdir, d=DATE, dom=DOMAINS, ksize=KSIZES),
        expand("{o}/lineages.{dom}.csv", o=outdir, dom=DOMAINS),
        expand("{o}/report/report.{d}-{dom}.html", o=outdir, d=DATE, dom=DOMAINS),

rule step_one:
    input:
        expand("{o}/lineages.{dom}.csv", o=outdir, dom=DOMAINS),


rule build_genbank:
    input:
        expand("{o}/genbank-{d}-{dom}-k{ksize}.zip", o=outdir, d=DATE, dom=DOMAINS, ksize=KSIZES),

rule clean_genbank:
    input:
        expand("{o}/genbank-{d}-{dom}-k{ksize}.clean.zip", o=outdir, d=DATE, dom=DOMAINS, ksize=KSIZES),

rule missing_genbank:
    input:
        expand("{o}/genbank-{d}-{dom}-k{ksize}.missing.zip", o=outdir, d=DATE, dom=DOMAINS, ksize=KSIZES),

rule check_genbank:
    input:
        expand("{o}/data/genbank-{d}-{dom}-k{ksize}.zip.check", o=outdir, d=DATE, dom=DOMAINS, ksize=KSIZES),

rule tax_genbank:
    input:
        expand("{o}/lineages.{dom}.csv", o=outdir, dom=DOMAINS),


#### Target Rules ####


rule download_assembly_summary:
    output:
        good = '{o}/data/assembly_summary.{dom}.txt',
    shell: """
        url_good="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.dom}/assembly_summary.txt"

        server_status_good=$(curl -L -o /dev/null -w "%{{http_code}}" -s "$url_good")

        if [ "$server_status_good" -eq 200 ]; then
            curl -L "$url_good" > {output.good}
        else
            echo "Failed to download files"
            echo "Server status code $server_status_good"
        fi
    """

rule download_historical_summary:
    output:
        bad = '{o}/data/assembly_summary_historical.{dom}.txt',
    shell: """
        url_bad="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.dom}/assembly_summary_historical.txt"

        server_status_bad=$(curl -L -o /dev/null -w "%{{http_code}}" -s "$url_bad")

        if [ "$server_status_bad" -eq 200 ]; then
            curl -L "$url_bad" > {output.bad}
        else
            echo "Failed to download files"
            echo "Server status code $server_status_bad"
        fi
    """

# create a 1% sub_assembly file for testing!!!
rule test_with_sub_assembly_summary:
    input:
        good = '{o}/data/assembly_summary.{dom}.txt',
        bad = '{o}/data/assembly_summary_historical.{dom}.txt',
    output:
        good = '{o}/data/sub_assembly_summary.{dom}.txt',
        bad = '{o}/data/sub_assembly_summary_historical.{dom}.txt',
    shell:"""
        echo {input.good}
        cat {input.good} | wc -l
        awksize 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01 || FNR<4) print $0}}' {input.good} > {output.good}
        cat {output.good} | wc -l

        echo {input.bad}
        cat {input.bad} | wc -l
        awksize 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01 || FNR<4) print $0}}' {input.bad} > {output.bad}
        cat {output.bad} | wc -l
    """

rule get_ss_db:
    params:
       old_date = OLD_DATES,
       indir = INDIR
    output:
        dbs = temporary(f"genbank-{OLD_DATES}-{{dom}}-k{{ksize}}.zip"),
    conda: "envs/sourmash.yaml"
    shell: """
            echo "Old date: {params.old_date}"

            echo "Checking if genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip exists..."
            if [ -e /group/ctbrowngrp/sourmash-db/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip ]; then

                echo "genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip exists!"
                echo "Linking existing file to $(pwd)"

                ln -s /group/ctbrowngrp/sourmash-db/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip $(pwd)/{output.dbs}
            else

                echo "File not found in primary directory. Checking backup directory {params.indir}..."

                if [ -e {params.indir}/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip ]; then
                    echo "File found in backup directory! Linking..."

                    ln -s {params.indir}/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip $(pwd)/{output.dbs}

                else

                    echo "File not found in backup directory. Attempting to download..."

                    curl -L https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-{params.old_date}/genbank-{params.old_date}-{wildcards.dom}-k{wildcards.ksize}.zip > {output.dbs}

                    if [ $? -ne 0 ]; then
                        echo "Download failed! Exiting with error."
                        exit 1
                    fi
                fi
            fi
   """

rule collect_all:
    input:
        unpack(createOldSingleManifest), #unpack(createSingleManifest)[0], #unpack the first of the keys i.e. original_dbs
    output:
        db = f"{{o}}/data/collect-mf.{OLD_DATES}-{{dom}}.csv",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 32 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig manifest --no-rebuild {input.original_dbs} -o {output.db}
    """

rule cleanse_manifest:
    input:
        unpack(getInputFilesForManifest),
        script = "scripts/update_sourmash_dbs.py",
        manifest = f"{{o}}/data/collect-mf.{OLD_DATES}-{{dom}}.csv",
    output:
        clean = "{o}/data/mf-clean.{d}-{dom}.csv",
        reversion = "{o}/data/updated-versions.{d}-{dom}.csv",
        report = "{o}/data/update-report.{d}-{dom}.txt",
        missing = "{o}/data/missing-genomes.{d}-{dom}.csv",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean} --updated-version {output.reversion} --report {output.report} --missing-genomes {output.missing}
    """

rule picklist_clean_db:
    benchmark: '{o}/benchmark/picklist_clean_db.{d}-{dom}-k{ksize}.tsv'
    input:
        clean = "{o}/data/mf-clean.{d}-{dom}.csv",
        dbs = f"genbank-{OLD_DATES}-{{dom}}-k{{ksize}}.zip",
    output:
        woohoo = temporary("{o}/genbank-{d}-{dom}-k{ksize}.clean.zip"),
    conda: "envs/sourmash.yaml",
    params:
        old_date = OLD_DATES
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=10,
        partition="bmh",
    shell:'''
        echo "Cleaning {input.dbs}..."
        sourmash sig cat {input.dbs} --picklist {input.clean}:name:name -ksize {wildcards.ksize} -o {output.woohoo}
        echo "{input.dbs} cleaned and stored as {output.woohoo}"
    '''

rule gather_sketch_reversioned:
    benchmark: '{o}/benchmark/gather_sketch_reversioned.{d}-{dom}.tsv'
    input:
        reversion = "{o}/data/updated-versions.{d}-{dom}.csv",
    output:
        failed = "{o}/data/update.{d}-{dom}.failures.csv",
        db = temporary("{o}/genbank-{d}-{dom}.rever.zip"),
    conda: "envs/directsketch.yaml"
    threads: 32
    resources:
        mem_mb = 100 * 1024,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=50,
        partition="bmh",
    params:
        k_list = lambda wildcards: ",".join([f"ksize={ksize}" for ksize in KSIZES]),
        #k_list = lambda wildcards: f"ksize={','.join([f'{ksize}' for ksize in KSIZES])}",
        scale = config.get('scale_value'),
        api_key = NCBI_API_KEY,
        threads = lambda wildcards: int(10 if NCBI_API_KEY and NCBI_API_KEY.strip() else 3)
    log:
        "logs/gather_sketch_reversioned.{o}_{d}_{dom}.log"
    shell:'''
        sourmash scripts gbsketch {input.reversion} -o {output.db} --failed {output.failed} \
            --param-string "dna,{params.k_list},scaled={params.scale},abund" \
            -a '{params.api_key}' -r 10 -n {params.threads} -g 2> {log}
    '''

rule cat_to_clean_reversioned:
    input:
        dir = "{o}/genbank-{d}-{dom}.rever.zip",
        db = "{o}/genbank-{d}-{dom}-k{ksize}.clean.zip",
        missing = "{o}/data/update.{d}-{dom}.failures.csv",
    output:
        woohoo = temporary("{o}/genbank-{d}-{dom}-k{ksize}.update.zip"),
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=100,
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig cat {input.dir} {input.db} -ksize {wildcards.ksize} -o {output.woohoo}
    """

rule gather_sketch_missing:
    benchmark: '{o}/benchmark/gather_sketch_missing.{d}-{dom}.tsv'
    input:
         missing = "{o}/data/missing-genomes.{d}-{dom}.csv",
    output:
        failed = "{o}/data/missing-genomes.{d}-{dom}.failures.csv",
        db = temporary("{o}/genbank-{d}-{dom}.miss.zip"),
    conda: "envs/directsketch.yaml"
    resources:
        mem_mb = 128 * 1024,
        time = lambda wildcards, attempt: 10 * 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 10 * 24 * 60 * attempt,
        allowed_jobs=50,
        partition="bmh",
    threads: 32
    params:
        k_list = lambda wildcards: ",".join([f"ksize={ksize}" for ksize in KSIZES]),
        #k_list = lambda wildcards: f"ksize={','.join([f'{ksize}' for ksize in KSIZES])}",
        scale = config.get('scale_value'),
        log = "logs/gather_sketch_missing.{d}_{dom}.log",
        api_key = NCBI_API_KEY,
        threads = lambda wildcards: int(10 if NCBI_API_KEY and NCBI_API_KEY.strip() else 3)
    shell:'''
        sourmash scripts gbsketch {input.missing} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled={params.scale},abund" \
            -a '{params.api_key}' -r 10 -n {params.threads} -g 2> {params.log}
    '''

rule cat_to_clean_missing:
    input:
        rever = "{o}/genbank-{d}-{dom}.rever.zip",
        miss = "{o}/genbank-{d}-{dom}.miss.zip",
        clean = "{o}/genbank-{d}-{dom}-k{ksize}.clean.zip",
    output:
        woohoo = "{o}/genbank-{d}-{dom}-k{ksize}.zip",
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
        time = lambda wildcards, attempt: 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 24 * 60 * attempt,
        allowed_jobs=100,
        partition="bmm",
    shell: """
        sourmash sig cat {input.miss} {input.rever} {input.clean} -ksize {wildcards.ksize} -o {output.woohoo}
    """

rule collect_complete:
    input:
        unpack(createNewSingleManifest), #unpack the second manifest i.e. new_dbs
    output:
        db = f"{{o}}/data/collect-mf.{{d}}-{{dom}}.csv",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 32 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig manifest --no-rebuild {input} -o {output.db}
    """

rule picklist_check:
    input:
        dbs_manifest = "{o}/data/collect-mf.{d}-{dom}.csv",
        tax_picklist = '{o}/lineages.{dom}.csv',
    output:
        missing = "{o}/data/genbank-{d}-{dom}.missing.csv",
        manifest = "{o}/data/genbank-{d}-{dom}.existing.csv",
    params:
        log = "logs/genbank-{d}-{dom}.picklist_check.log",
        first_ksize = lambda wildcards: KSIZES[0],
    conda: "envs/sourmash.yaml"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: 6 * 1024 * attempt,
        time= 10000,
        partition='high2',
    shell:
        """
        sourmash sig checksize -ksize {params.first_ksize} \
            --picklist {input.tax_picklist}:ident:ident \
            {input.dbs_manifest} --output-missing {output.missing} \
            --save-manifest {output.manifest} 2> {params.log}
        touch {output.missing}
        """

rule make_manual_files:
    input:
        script = "scripts/gather_failed.sh",
        missing = "{o}/data/missing-genomes.{d}-{dom}.csv",
        reversion = "{o}/data/updated-versions.{d}-{dom}.csv",
        checksize = "{o}/data/genbank-{d}-{dom}.missing.csv",
    output:
        output = "{o}/workflow-cleanup/manual-download.{d}-{dom}.csv",
        manual = "{o}/workflow-cleanup/manual-check.{d}-{dom}.csv",
        log = "{o}/workflow-cleanup/log.{d}-{dom}.txt",
    shell:"""
        {input.script} {wildcards.o} {wildcards.d} {wildcards.dom} 2>&1 | tee {output.log}
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
        "{o}/data/assembly_summary.{dom}.txt",
        "taxdump/nodes.dmp",
        "taxdump/names.dmp",
        "scripts/make-lineage-csv.py",
        "scripts/ncbi_taxdump_utils.py",
    output:
        "{o}/lineages.{dom}.csv"
    params:
        ictv_cmd = lambda w: " --ictv " if 'viral' in w.dom else '',
    shell:
        "python scripts/make-lineage-csv.py taxdump/{{nodes.dmp,names.dmp}} {input[0]} -o {output} {params.ictv_cmd}"

### create a report with sourmash sig summarize for the databases... and sourmash compare(?)

rule quarto_report:
    input:
        unpack(getInputFilesForManifest),
        report = "{o}/data/update-report.{d}-{dom}.txt",
        #new_mf = f"{{o}}/data/collect-mf.{d}-{{dom}}.csv",
        #new_mf = f"{{o}}/data/collect-mf.{str(DATE)}-{{dom}}.csv",
        new_mf = "{o}/data/collect-mf.{d}-{dom}.csv",
        old_mf = f"{{o}}/data/collect-mf.{OLD_DATES}-{{dom}}.csv",
        failures = "{o}/data/missing-genomes.{d}-{dom}.failures.csv",
        missing = "{o}/data/genbank-{d}-{dom}.missing.csv",
        gathered = "{o}/data/genbank-{d}-{dom}.existing.csv",
        lineage = "{o}/lineages.{dom}.csv",
        recovered = "{o}/workflow-cleanup/log.{d}-{dom}.txt",
    output:
        "{o}/report/report.{d}-{dom}.html",
    params:
        log = "logs/{d}_{dom}_report.log",
        report_title = "Genbank's {dom} Database Update Report",
        old_date = OLD_DATES,
        old_db = lambda wildcards: ",".join([f'"genbank-{OLD_DATES}-{wildcards.dom}-k{ksize}.zip"' for ksize in KSIZES]),
        new_db = lambda wildcards: ",".join([f'"{wildcards.o}/genbank-{wildcards.d}-{wildcards.dom}-k{ksize}.zip"' for ksize in KSIZES]),
        man = "{o}/workflow-cleanup/manual-download.{d}-{dom}.csv",
        man_check = "{o}/workflow-cleanup/manual-check.{d}-{dom}.csv",
        man_out = "{o}/workflow-cleanup/manual-download.{d}-{dom}.zip",
        man_fail = "{o}/workflow-cleanup/manual-download.{d}-{dom}.failed.csv",
        man_log = "{o}/workflow-cleanup/manual-download.{d}-{dom}.log",
        k_list = lambda wildcards: ",".join([f"ksize={ksize}" for ksize in KSIZES]),
        scale = config.get('scale_value'),
    conda: "envs/quarto.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell:
        """
        # Mimicing https://github.com/ETH-NEXUS/quarto_example/blob/main/workflow/rules/clean_data_report.smk
        # This will be a stand-alone html document and needs to embed-resources
        # https://quarto.org/docs/output-formats/html-publishing.html#standalone-html
        mkdir -p {wildcards.d}-{wildcards.dom}.temp
        cp scripts/report.qmd {wildcards.d}-{wildcards.dom}.temp/
        cd {wildcards.d}-{wildcards.dom}.temp/

        DIRNAME=$(dirname "{output}")

        quarto render report.qmd \
            -P details_files:{input.report} -P old_details:{params.old_date} \
            -P new_details:{wildcards.d} -P old_db:{params.old_db} \
            -P new_db:{params.new_db} -P old_mf:{input.old_mf} \
            -P new_mf:{input.new_mf} -P failures:{input.failures} \
            -P report_title:"{params.report_title}" -P config:"{config}" \
            -P good_assm:{input.good} -P bad_assm:{input.bad} \
            -P missed:{input.missing} -P gathered:{input.gathered} \
            -P lineage:{input.lineage} -P recovered:{input.recovered} \
            -P manual_download:{params.man} -P manual_output:{params.man_out} \
            -P manual_failed:{params.man_fail} -P manual_log:{params.man_log} \
            -P manual_check:{params.man_check} -P output_dir:{wildcards.o} \
            -P k_list:{params.k_list} -P scale:{params.scale}
        # 2> {params.log}

        # the `--output` arg adds an unnecessary `../` to the output file path
        # https://github.com/quarto-dev/quarto-cli/issues/10129
        # just mving it back in place

        mv report.html {output}
        cd ..
        rm -rf {wildcards.d}-{wildcards.dom}.temp
        """
