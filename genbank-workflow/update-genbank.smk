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

# Dictionary for dynamic slurm batch allocations with correct resources
#PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 33], 4: ['med2', 33], 5: ['high2', 100]}
PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 33], 4: ['bmm', 33], 5: ['bmh', 100]}

# Create a file dictionary for normal and test runs for rule cheat_mainfest
def getInputFilesForManifest(wildcards):
    files = dict()
    if config.get('output_directory') == 'test':
        files["good"] = f"{wildcards.o}/data/sub_assembly_summary.{wildcards.D}.txt"
        files["bad"] = f"{wildcards.o}/data/sub_assembly_summary_historical.{wildcards.D}.txt"
    else:
        files["good"] = f"{wildcards.o}/data/assembly_summary.{wildcards.D}.txt"
        files["bad"] = f"{wildcards.o}/data/assembly_summary_historical.{wildcards.D}.txt"
    return files

rule all:
    input:
        expand("{o}/genbank-{d}-{D}-k{k}.zip", o=outdir, d=DATE, D=DOMAINS, k=KSIZES),
        expand("{o}/lineages.{D}.csv", o=outdir, D=DOMAINS),
        expand("{o}/report/report.{d}-{D}.html", o=outdir, d=DATE, D=DOMAINS),

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
    shell: """
        url_good="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary.txt"

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
        bad = '{o}/data/assembly_summary_historical.{D}.txt',
    shell: """
        url_bad="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.D}/assembly_summary_historical.txt"

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
        good = '{o}/data/assembly_summary.{D}.txt',
        bad = '{o}/data/assembly_summary_historical.{D}.txt',
    output:
        good = '{o}/data/sub_assembly_summary.{D}.txt',
        bad = '{o}/data/sub_assembly_summary_historical.{D}.txt',
    shell:"""
        echo {input.good}
        cat {input.good} | wc -l
        awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01 || FNR<4) print $0}}' {input.good} > {output.good}
        cat {output.good} | wc -l

        echo {input.bad}
        cat {input.bad} | wc -l
        awk 'BEGIN {{srand()}} !/^$/ {{ if (rand() <= .01 || FNR<4) print $0}}' {input.bad} > {output.bad}
        cat {output.bad} | wc -l
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

rule collect_all:
    input:
        dbs = lambda wildcards: expand(f"genbank-{OLD_DATES}-{{D}}-k{{k}}.zip", D = [wildcards.D], k = KSIZES[0])
    output:
        db = f"{{o}}/data/collect-mf.{OLD_DATES}-{{D}}.csv",
    params:
        first_k = lambda wildcards: KSIZES[0],
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 32 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig manifest --no-rebuild {input.dbs} -o {output.db}
    """

rule cleanse_manifest:
    input:
        unpack(getInputFilesForManifest),
        script = "scripts/update_sourmash_dbs.py",
        manifest = f"{{o}}/data/collect-mf.{OLD_DATES}-{{D}}.csv",
    output:
        clean = "{o}/data/mf-clean.{d}-{D}.csv",
        reversion = "{o}/data/updated-versions.{d}-{D}.csv",
        report = "{o}/data/update-report.{d}-{D}.txt",
        missing = "{o}/data/missing-genomes.{d}-{D}.csv",
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
    input:
        clean = "{o}/data/mf-clean.{d}-{D}.csv",
        dbs = f"genbank-{OLD_DATES}-{{D}}-k{{k}}.zip",
    output:
        woohoo = temporary("{o}/genbank-{d}-{D}-k{k}.clean.zip"),
    conda: "envs/sourmash.yaml",
    params:
        old_date = OLD_DATES
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=10,
        partition="bmh",
    benchmark:
        'benchmarks/picklist_clean_db.{o}-{d}-{D}-k{k}.tsv'
    shell:'''
        echo "Cleaning {input.dbs}..."
        sourmash sig cat {input.dbs} --picklist {input.clean}:name:name -k {wildcards.k} -o {output.woohoo}
        echo "{input.dbs} cleaned and stored as {output.woohoo}"
    '''

rule gather_sketch_reversioned:
    input:
        reversion = "{o}/data/updated-versions.{d}-{D}.csv",
    output:
        failed = "{o}/data/update.{d}-{D}.failures.csv",
        db = temporary("{o}/genbank-{d}-{D}.rever.zip"),
    conda: "envs/directsketch.yaml"
    benchmark:
        'benchmarks/gather_sketch_reversioned.{o}-{d}-{D}.tsv'
    threads: 3
    resources:
        mem_mb = 100 * 1024,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=100,
        partition="bmh",
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
    log:
        "logs/gather_sketch_reversioned.{o}_{d}_{D}.log"
    shell:'''
        sourmash scripts gbsketch {input.reversion} -o {output.db} --failed {output.failed} \
            --param-string "dna,{params.k_list},scaled=1000,abund" -r 5 -g 2> {log}
    '''

rule cat_to_clean_reversioned:
    input:
        dir = "{o}/genbank-{d}-{D}.rever.zip",
        db = "{o}/genbank-{d}-{D}-k{k}.clean.zip",
        missing = "{o}/data/update.{d}-{D}.failures.csv",
    output:
        woohoo = temporary("{o}/genbank-{d}-{D}-k{k}.update.zip"),
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=100,
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    benchmark:
        'benchmarks/cat_to_clean_reversioned.{o}-{d}-{D}-k{k}.tsv'
    shell: """
        sourmash sig cat {input.dir} {input.db} -k {wildcards.k} -o {output.woohoo}
    """

rule gather_sketch_missing:
    input:
         missing = "{o}/data/missing-genomes.{d}-{D}.csv",
    output:
        failed = "{o}/data/missing-genomes.{d}-{D}.failures.csv",
        db = temporary("{o}/genbank-{d}-{D}.miss.zip"),
    conda: "envs/directsketch.yaml"
    resources:
        mem_mb = 100 * 1024,
        time = lambda wildcards, attempt: 168 * 60 * attempt,
        runtime = lambda wildcards, attempt: 168 * 60 * attempt,
        allowed_jobs=100,
        partition="bmh",
    threads: 3
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        #k_list = lambda wildcards: f"k={','.join([f'{k}' for k in KSIZES])}",
        log = "logs/gather_sketch_missing.{d}_{D}.log"
    shell:'''
        sourmash scripts gbsketch {input.missing} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled=1000,abund" -r 5 -g 2> {params.log}
    '''

# Could this be changed to `sourmash sig check` instead?
checkpoint cat_to_clean_missing:
    input:
        rever = "{o}/genbank-{d}-{D}.rever.zip",
        miss = "{o}/genbank-{d}-{D}.miss.zip",
        db = "{o}/genbank-{d}-{D}-k{k}.clean.zip",
    output:
        woohoo = protected("{o}/genbank-{d}-{D}-k{k}.zip"),
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=100,
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig cat {input.miss} {input.rever} {input.db} -k {wildcards.k} -o {output.woohoo}
    """

rule collect_complete:
    input:
        dbs = lambda wildcards: expand("{o}/genbank-{d}-{D}-k{k}.zip", o = [wildcards.o], d = [wildcards.d], D = [wildcards.D], k = KSIZES[0])
    output:
        db = "{o}/data/collect-mf.{d}-{D}.csv",
    params:
        first_k = lambda wildcards: KSIZES[0],
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 32 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig manifest --no-rebuild {input.dbs} -o {output.db}
    """

rule picklist_check:
    input:
        dbs_manifest = "{o}/data/collect-mf.{d}-{D}.csv",
        tax_picklist = '{o}/lineages.{D}.csv',
    output:
        missing = "{o}/data/genbank-{d}-{D}.missing.csv",
        manifest = "{o}/data/genbank-{d}-{D}.existing.csv",
    params:
        log = "logs/genbank-{d}-{D}.picklist_check.log",
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

### create a report with sourmash sig summarize for the databases... and sourmash compare(?)

rule quarto_report:
    input:
        unpack(getInputFilesForManifest),
        report = "{o}/data/update-report.{d}-{D}.txt",
        new_mf = "{o}/data/collect-mf.{d}-{D}.csv",
        old_mf = f"{{o}}/data/collect-mf.{OLD_DATES}-{{D}}.csv",
        failures = "{o}/data/missing-genomes.{d}-{D}.failures.csv",
        missing = "{o}/data/genbank-{d}-{D}.missing.csv",
        gathered = "{o}/data/genbank-{d}-{D}.existing.csv",
        lineage = "{o}/lineages.{D}.csv"
    output:
        "{o}/report/report.{d}-{D}.html"
    params:
        log = "../logs/{d}_{D}_report.log",
        report_title = "Genbank's {D} Database Update Report",
        old_date = OLD_DATES,
        old_db = lambda wildcards: ",".join([f'"genbank-{OLD_DATES}-{wildcards.D}-k{k}.zip"' for k in KSIZES]),
        new_db = lambda wildcards: ",".join([f'"{wildcards.o}/genbank-{wildcards.d}-{wildcards.D}-k{k}.zip"' for k in KSIZES]),
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
        mkdir -p {wildcards.d}-{wildcards.D}.temp
        cp scripts/report.qmd {wildcards.d}-{wildcards.D}.temp/
        cd {wildcards.d}-{wildcards.D}.temp/

        DIRNAME=$(dirname "{output}")

        quarto render report.qmd \
            -P details_files:{input.report} -P old_details:{params.old_date} \
            -P new_details:{wildcards.d} -P old_db:{params.old_db} \
            -P new_db:{params.new_db} -P old_mf:{input.old_mf} \
            -P new_mf:{input.new_mf} -P failures:{input.failures} \
            -P report_title:"{params.report_title}" -P config:"{config}" \
            -P good_assm:{input.good} -P bad_assm:{input.bad} \
            -P missed:{input.missing} -P gathered:{input.gathered} \
            -P lineage:{input.lineage}
        # 2> {params.log}

        # the `--output` arg adds an unnecessary `../` to the output file path
        # https://github.com/quarto-dev/quarto-cli/issues/10129
        # just mving it back in place

        mv report.html {output}
        cd ..
        rm -rf {wildcards.d}-{wildcards.D}.temp
        """
