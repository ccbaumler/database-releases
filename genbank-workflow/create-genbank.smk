###
# This workflow will create a genbank database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s create-genbank.smk -j 6 --use-conda --retries=3 --rerun-incomplete --resources allowed_jobs=100 --latency-wait=60
# Jobs (`-j`) should be 3 times the number of domains.
###
##

configfile: "config/create-genbank.yaml"

DATE = config.get("date")

if not DATE:
    import time
    DATE = [time.strftime("%Y%m%d")]

outdir = [
    f"{config.get('output_directory') if config.get('output_directory') is not None else '..'}/genbank-{date}"
    for date in DATE]

NEW_DOMAINS = config.get('domains')

KSIZES = config.get('k_values')

email = config.get('email')

if email:
    onsuccess:
        print("\nWorkflow finished without error\n")
        shell("mail -s 'Workflow finished without error' {email} < {log}")
    
    onerror:
        print("\nAn error occurred\n")
        shell("mail -s 'an error occurred' {email} < {log}")

wildcard_constraints:
    k = "\d{2}+",
    ND = "\w[^-]+",
    d = "\d+",

# Create a file dictionary for normal and test runs for rule cheat_mainfest
def getInputFilesForManifest(wildcards):
    files = dict()
    if config.get('output_directory') == 'test':
        files["good"] = f"{wildcards.o}/data/sub_assembly_summary.{wildcards.ND}.txt"
        files["bad"] = f"{wildcards.o}/data/sub_assembly_summary_historical.{wildcards.ND}.txt"
    else:
        files["good"] = f"{wildcards.o}/data/assembly_summary.{wildcards.ND}.txt"
        files["bad"] = f"{wildcards.o}/data/assembly_summary_historical.{wildcards.ND}.txt"
    return files

rule all:
    input:
        expand("{o}/genbank-{d}-{ND}-k{k}.zip", o=outdir, d=DATE, ND=NEW_DOMAINS, k=KSIZES),
        expand("{o}/lineages.{ND}.csv", o=outdir, ND=NEW_DOMAINS),

rule build_databases:
    input:
        expand("{o}/genbank-{d}-{ND}-k{k}.zip", o=outdir, d=DATE, ND=NEW_DOMAINS, k=KSIZES),

rule build_lineages:
    input:
        expand("{o}/lineages.{ND}.csv", o=outdir, ND=NEW_DOMAINS),

rule download_assembly_summary:
    output:
        good = '{o}/data/assembly_summary.{ND}.txt',
        bad = '{o}/data/assembly_summary_historical.{ND}.txt',
    shell: """
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.ND}/assembly_summary.txt > {output.good}
        curl -L https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{wildcards.ND}/assembly_summary_historical.txt > {output.bad} 
    """

# create a 1% sub_assembly file for testing!!!
rule test_with_sub_assembly_summary:
    input:
        good = '{o}/data/assembly_summary.{ND}.txt',
        bad = '{o}/data/assembly_summary_historical.{ND}.txt',
    output:
        good = '{o}/data/sub_assembly_summary.{ND}.txt',
        bad = '{o}/data/sub_assembly_summary_historical.{ND}.txt',
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

rule new_manifest:
    output:
        mani = '{o}/data/mf.{d}-{ND}-k{k}.csv',
    shell:"""
        cat > {output.mani} << EOF
# SOURMASH-MANIFEST-VERSION: 1.0
internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename
EOF
    """

rule cheat_manifest:
    input:
        unpack(getInputFilesForManifest),
        script = "scripts/update_sourmash_dbs.py",
        manifest = "{o}/data/mf.{d}-{ND}-k{k}.csv",
    output:
        clean = "{o}/data/mf-clean.{d}-{ND}-k{k}.csv",
        missing = "{o}/data/missing-genomes.{d}-{ND}-k{k}.csv",
    conda: "envs/sourmash.yaml",
    shell: """
        {input.script} {input.manifest} -a {input.good} -b {input.bad} -o {output.clean}  --missing-genomes {output.missing}
    """

rule check_txt_new:
    input:
        new = expand("{o}/data/missing-genomes.{d}-{ND}-k{k}.csv", o=outdir, d=DATE, ND=NEW_DOMAINS, k=KSIZES),
        script = "scripts/check_txt_files.py",
    output:
        solo = "{o}/data/missing.{d}-{ND,\w[^-]+}.csv",
    wildcard_constraints:
        ND="\w[^-]+"
    shell: """
        files=""
        for file in "{wildcards.o}/data/missing-genomes.*{wildcards.d}-{wildcards.ND}-k*.csv"; do
            files+=" $file"
        done
        {input.script} $files -o {output.solo}
    """

rule gather_sketch_db:
    input:
        solo = "{o}/data/missing.{d}-{ND}.csv",
    output:
        failed = "{o}/data/missing.{d}-{ND,\w[^-]+}.failures.csv",
        db = temporary("{o}/temp/genbank-{d}-{ND,\w[^-]+}.zip"),
    wildcard_constraints:
        ND="\w[^-]+"
    conda: "envs/directsketch.yaml",
    resources:
        mem_mb =  100 * 1024,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    threads: 3
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        scale = config.get('scale_value'),
        log = 'logs/gather_sketch_db.{d}_{ND}.log'
    shell:'''
        sourmash scripts gbsketch {input.solo} -o {output.db} --failed {output.failed} \
            --param-str "dna,{params.k_list},scaled={params.scale},abund" -r 1 2> {params.log}
    '''

rule extract_db:
    input:
        db = "{o}/temp/genbank-{d}-{ND}.zip",
    output:
        db = "{o}/genbank-{d}-{ND}-k{k}.zip",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = 80 * 1024,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    shell:"""
        sourmash signature cat {input.db} -k {wildcards.k} -o {output.db}
    """

# Include a quarto report showing all the info for the workflow

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
        "{o}/data/assembly_summary.{ND}.txt",
        "taxdump/nodes.dmp",
        "taxdump/names.dmp",
        "scripts/make-lineage-csv.py",
        "scripts/ncbi_taxdump_utils.py",
    output:
        "{o}/lineages.{ND}.csv"
    params:
        ictv_cmd = lambda w: " --ictv " if 'viral' in w.ND else '',
    shell:
        "python scripts/make-lineage-csv.py taxdump/{{nodes.dmp,names.dmp}} {input[0]} -o {output} {params.ictv_cmd}"
