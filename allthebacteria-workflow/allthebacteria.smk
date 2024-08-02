###
# This workflow will create an allthebacteria database for sourmash from the HPC cluster at UCDavis.
#
# To run:
# snakemake -s allthebacteria.smk -j 1 --use-conda --resources allowed_jobs=100
###

import os
import re
import random
import sys
import csv


configfile: "config/allthebacteria.yaml"

EMAIL=config.get('email')
RELEASES=config.get('releases')
FILES=config.get('files_of_paths')
KSIZES=config.get('k_values')
SCALE=config.get('scale_value')
OUTDIR = [
    f"{config.get('output_directory') if config.get('output_directory') is not None else '..'}/allthebacteria-r{release}"
    for release in RELEASES
    ]

wildcard_constraints:
    r = "\d+\.\d+", 

# Dictionary for dynamic slurm batch allocations with correct resources
#PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 33], 4: ['med2', 33], 5: ['high2', 100]}
PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 33], 4: ['bmm', 33], 5: ['bmh', 100]}

dir_name_dict = {}
pattern = r'r(\d+\.\d+)'
for file in FILES:
    match = re.search(pattern, file)
    if match:
        release = float(match.group(1))
        with open(file, 'r') as fp:
            file_list = fp.readlines()

        dir_name_list = [os.path.basename(f).split(".")[0] for f in file_list]

        if config.get('output_directory') == 'test':
            selected_lines = file_list[-2:] #random.sample(file_list, 2) #random.choice(file_list)
            dir_name_list = [os.path.basename(line).split(".")[0] for line in selected_lines]

            directory_path = f"test/allthebacteria-r{release}"

            if not os.path.exists(directory_path):
                os.makedirs(directory_path)
            else:
                print(f"The directory '{directory_path}' already exists.")
    
            test_filename = f"{directory_path}/sub-{file}"
            
            with open(test_filename, 'w') as fp:
                fp.write(''.join(selected_lines))

            print(f"Selected line for test: {dir_name_list} written to {test_filename}")

        else:

            if not os.path.isdir(f'{OUTDIR[0]}'):
                os.makedirs(f'{OUTDIR[0]}')
            else:
                print(f"The directory '{OUTDIR[0]}' already exists.")

            filename = f"{OUTDIR[0]}/{file}"

            with open(filename, 'w') as fp:
                fp.write(''.join(file_list))

            print(f"File {file} copied to {OUTDIR[0]}")

        dir_name_dict[release] = dir_name_list
    else:
        raise ValueError(f"No release found in filename: {file}")

if EMAIL:
    onsuccess:
        print("\nWorkflow finished without error\n")
        shell("mail -s 'Workflow finished without error' {EMAIL} < {log}")
    
    onerror:
        print("\nAn error occurred\n")
        shell("mail -s 'an error occurred' {EMAIL} < {log}")

# Create a file dictionary for normal and test runs for rule download_links
def getInputFilesForLinks(wildcards):
    files = dict()
    for file in FILES:
        if config.get('output_directory') == 'test':
            files["data"] = f"{wildcards.o}/sub-{file}"
        else:
            files["data"] = f"{wildcards.o}/{file}"
        return files

rule all:
    input:
        expand("{o}/allthebacteria-r{r}-sigs/{dir_name}/{dir_name}.zip", o=OUTDIR, r=RELEASES, dir_name=[dn for r in RELEASES for dn in dir_name_dict[r]]),
        expand("{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}-k{k}.zip", o=OUTDIR, r=RELEASES, k=KSIZES),
        expand("{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-mf.csv", o=OUTDIR, r=RELEASES),
        expand("{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-missing.csv", o=OUTDIR, r=RELEASES),
        expand("{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-existing.csv", o=OUTDIR, r=RELEASES),
        expand("{o}/report/report.{r}.html", o=OUTDIR, r=RELEASES),

rule check_sketch:
    input:
        zip_file = expand("{o}/allthebacteria-r{r}-sigs/{dir_name}/{dir_name}.zip", o=OUTDIR, r=RELEASES, dir_name=[dn for r in RELEASES for dn in dir_name_dict[r]]),
    log:
        '{o}/logs/check_sketch/allthebacteria-r{r}-sig-{dir_name}.log'
    shell:"""
        unzip -v {input.zip_file} 2>&1 | tee -a {log}
    """

rule report:
    input:
        expand("{o}/report/report.allthebacteria-r{r}.html", o=OUTDIR, r=RELEASES),

rule make_links:
    input:
        script = "scripts/ftp_link_list.py",
    output:
        meta = "{o}/allthebacteria-r{r}-metadata-links.txt",
        data = "{o}/allthebacteria-r{r}-data-links.txt",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
	partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        {input.script} ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/{wildcards.r}/metadata -s ftp -p {output.meta}

        {input.script} ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/{wildcards.r}/assembly -s ftp -p {output.data}
    """

rule download_links:
    input:
        unpack(getInputFilesForLinks),
        meta = "{o}/allthebacteria-r{r}-metadata-links.txt",
    output:
        tsv = "{o}/allthebacteria-r{r}-metadata/sylph.tsv.gz",
        txt = "{o}/allthebacteria-r{r}-metadata/sample_list.txt.gz",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
	partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    threads: 
        10
    shell: """
        # read the txt files, use only 1 link (-n), parallelize wget cmd upto 10 times (-P), path for downloads (-P), if interrupted continue the download (--continue), do not download anything that exists (--no-clobber), try up to 3 times per download (--tries), wait 1 sec per try (--wait), quiet output (-nv) 
        cat {input.data} | xargs -n 1 -P {threads} wget -P {wildcards.o}/allthebacteria-r{wildcards.r}-data/ --continue --no-clobber --tries=3 --wait=1 -nv && cat {input.meta} | xargs -n 1 -P {threads} wget -P {wildcards.o}/allthebacteria-r{wildcards.r}-metadata/ --continue --no-clobber --tries=3 --wait=1 -nv
    """

rule extract_info:
    input:  
        info = "{o}/allthebacteria-r{r}-metadata/sylph.tsv.gz",
        samples = "{o}/allthebacteria-r{r}-metadata/sample_list.txt.gz",
    output:
        info = "{o}/allthebacteria-r{r}-metadata/sylph.tsv",
        samples = "{o}/allthebacteria-r{r}-metadata/sample_list.txt",
        idents = "{o}/allthebacteria-r{r}-metadata/ident_list.txt",
    resources:
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell:"""
        gzip -kvdf {input.info}
        gzip -kvdf {input.samples}

        echo -e "ident" | cat - {output.samples} > {output.idents}
    """

rule extract_files:
    input:
        unpack(getInputFilesForLinks),
        placeholder = "{o}/allthebacteria-r{r}-metadata/sylph.tsv",
        script = "scripts/extract_check.sh",
    output:
        seq_dir = temporary(directory("{o}/allthebacteria-r{r}-seqs/{dir_name}")),
    conda: "envs/branchwater.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs = 10, #lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
#    group: "group"
    params:
        outdir = "{o}/allthebacteria-r{r}-seqs",
        data = "{o}/allthebacteria-r{r}-data/{dir_name}.asm.tar.xz", # needs to be in params because smk is looking for the file in upstream rule output
        download_dir = "{o}/allthebacteria-r{r}-data",
    shell:
        """
        echo "Extracting sequence files from compressed archives"

        # a script to double-check the failed extracts and re-download/extract if they were corrupted
        {input.script} -f {params.data} -o {params.outdir} -d {wildcards.dir_name} -p {params.download_dir} -i {input.data}

        #pv {params.data} | tar --use-compress-program="xz -T0 -q" --skip-old-files -xf - -C {params.outdir}
        #tar --skip-old-files -xvf {params.data} -C {params.outdir}
        echo "Extracted $(ls -1 {output.seq_dir} | wc -l) files!"
        """

rule build_csv:
    input:
        seq_dir = "{o}/allthebacteria-r{r}-seqs/{dir_name}",
        info_file = "{o}/allthebacteria-r{r}-metadata/sylph.tsv",
    output:
        csv_file = "{o}/allthebacteria-r{r}-sigs/{dir_name}/manysketch.csv",
        missing_file = touch("{o}/allthebacteria-r{r}-sigs/{dir_name}/missing.csv"),
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    params:
        dir_name = "{dir_name}"
    run:
        def load_info_dict(info_files=[]):
            info_dict = {}
            for info_file in info_files:
                with open(info_file, 'r') as fp:
                    reader = csv.DictReader(fp, delimiter='\t')
                    for row in reader:
                        info_dict[row['sample']] = row['sample'] + " " + row['Contig_name']
            return info_dict

        info_dict = load_info_dict([f"{wildcards.o}/allthebacteria-r{release}-metadata/sylph.tsv" for release in RELEASES],)

        print(f"Processing {input.seq_dir}")

        with open(output.csv_file, 'w', newline='') as csvfile, open(output.missing_file, 'w') as missing_fp:
             csvwriter = csv.writer(csvfile)
             csvwriter.writerow(['name', 'genome_filename', 'protein_filename'])

             filepaths = [os.path.join(input.seq_dir, f) for f in os.listdir(input.seq_dir) if os.path.isfile(os.path.join(input.seq_dir, f))]
             num_files = len(filepaths) - 1

             increment = int(num_files * 0.05)
             target_line = increment

             line_count = 0

             for filepath in filepaths:
                 filename_ext = os.path.basename(filepath)
                 filename = os.path.splitext(filename_ext)[0]

                 if filename in info_dict:
                     name = info_dict[filename]
                     name = name.replace(',', '')

                     csvwriter.writerow([name, filepath,''])

                     line_count += 1

                     if line_count >= target_line:
                         print(f"Progress: {line_count}/{num_files} lines written")
                         target_line += increment
                 else:
                     if not filename.startswith('.'):
                         missing_fp.write(f"{filename, filepath}\n")

                         csvwriter.writerow([filename, filepath,''])

                     line_count += 1

                     if line_count >= target_line:
                         print(f"Progress: {line_count}/{num_files} lines written")
                         target_line += increment

rule sketch_seqs:
    input:
        csv_file = "{o}/allthebacteria-r{r}-sigs/{dir_name}/manysketch.csv",
        seq_dir = "{o}/allthebacteria-r{r}-seqs/{dir_name}",
    output:
        zip_file = temporary("{o}/allthebacteria-r{r}-sigs/{dir_name}/{dir_name}.zip"),
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    threads: 32
    conda: 'envs/branchwater.yaml'
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        scale = {SCALE},
    shell:"""
        sourmash scripts manysketch {input.csv_file} -p {params.k_list},scaled={params.scale},abund -o {output.zip_file}

        # To prevent 'sourmash.exceptions.Panic: sourmash panicked: thread 'unnamed' panicked with 'called `Result::unwrap()` on an `Err` value: InvalidArchive("Couldn't find End Of Central Directory Record")' at src/core/src/storage.rs:367'
        UNZIP_OUTPUT=$(unzip -v {output.zip_file})
        UNZIP_EXIT_CODE=$?

        if [ $UNZIP_EXIT_CODE -ne 0 ]; then
            echo "Error found in sketch file. Checking the specific error..."

            if grep -q "End-of-central-directory" "$UNZIP_OUTPUT"; then
                echo "Detected 'End-of-central-directory' error. Removing corrupted sketch and re-sketching..."

                rm {output.zip_file}
                sourmash scripts manysketch {input.csv_file} -p {params.k_list},scaled={params.scale},abund -o {output.zip_file}

                UNZIP_OUTPUT=$(unzip -v {output.zip_file})
                UNZIP_EXIT_CODE=$?

                if [ $UNZIP_EXIT_CODE -ne 0 ]; then
                    echo "Sketching {output.zip_file} failed again."
                    echo "Please check extracted sequences in {input.seq_dir}"
                    echo "Or re-extract {wildcards.o}/allthebacteria-r{wildcards.r}-data/{wildcards.dir_name}.asm.tar.xz"
                    exit 1
                else
                    echo "Sketching completed successfully after second attempt."
                fi
            else
                echo "Sketching failed due to different error. Exiting."
                exit 1
            fi
        else
            echo "Sketching completed correctly."
        fi
    """

rule cat_all:
    input:
        zip_file = expand("{o}/allthebacteria-r{r}-sigs/{dir_name}/{dir_name}.zip", o=OUTDIR, r=RELEASES, dir_name=[dn for r in RELEASES for dn in dir_name_dict[r]]),
    output:
        all = "{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}.zip",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 5 * 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 5 * 24 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    conda: 'envs/sourmash.yaml',
    params:
        path = "{o}/allthebacteria-r{r}-sigs"
    shell:"""
        echo -e "\nFinding all sourmash signatures..."

        temp_file={params.path}/sig_list.txt
        find {params.path} -type f -name "*.zip" > $temp_file

        echo -e "\nFound $(wc -l $temp_file) signature files"
        echo -e "Combining all signature files now..."

        sourmash sig cat $(cat $temp_file) -o {output.all}
    """

rule manifest_manifest:
    input:
        db_zip = "{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}.zip",
    output:
        manifest = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-mf.csv",
    conda: 'envs/sourmash.yaml',
    shell:"""
        sourmash sig manifest --no-rebuild {input.db_zip} -o {output.manifest} 
    """

rule picklist_check:
    input:
        dbs_manifest = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-mf.csv",
        sample_picklist = '{o}/allthebacteria-r{r}-metadata/ident_list.txt',
    output:
        missing = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-missing.csv",
        manifest = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-existing.csv",
    params:
        log = "logs/allthebacteria-r{r}.picklist_check.log",
        first_k = KSIZES[0]
    conda: "envs/sourmash.yaml"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: 6 * 1024 * attempt,
        time= 10000,
        partition='high2',
    shell:
        """
        sourmash sig check -k {params.first_k} \
            --picklist {input.sample_picklist}:ident:ident \
            {input.dbs_manifest} --output-missing {output.missing} \
            --save-manifest {output.manifest} 2> {params.log}
        touch {output.missing}
        """

rule find_missing:
    input:
        script = "scripts/find_missing_files.sh",
        data_dir = "{o}/allthebacteria-r{r}-data/",
        missing = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-missing.csv",
        sample_picklist = '{o}/allthebacteria-r{r}-metadata/ident_list.txt',
    output:
        missing_files = '{o}/workflow-cleanup/allthebacteria-r{r}-missing-files.csv',
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    threads: 32
    shell:"""
        {input.script} {input.data_dir} {input.missing} {output.missing_files} {threads}
    """

rule cat_by_k:
    input:
        zip_file = "{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}.zip",
    output:
        k_zip_file = "{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}-k{k}.zip",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 5 * 24 * 60 * attempt,
        runtime = lambda wildcards, attempt: 5 *24 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    conda: 'envs/sourmash.yaml',
    shell:"""
        sourmash sig cat {input.zip_file} -k {wildcards.k} -o {output.k_zip_file}
    """

### create a report with sourmash sig summarize for the databases... and sourmash compare(?)

rule quarto_report:
    input:
        unpack(getInputFilesForLinks),
        meta = "{o}/allthebacteria-r{r}-metadata-links.txt",
        new_mf = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-mf.csv",
        names = "{o}/allthebacteria-r{r}-metadata/sylph.tsv",
        idents = "{o}/allthebacteria-r{r}-metadata/ident_list.txt",
        missing = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-missing.csv",
        existing = "{o}/allthebacteria-r{r}-metadata/allthebacteria-r{r}-existing.csv",
        missing_files = "{o}/workflow-cleanup/allthebacteria-r{r}-missing-files.csv",
    output:
        "{o}/report/report.allthebacteria-r{r}.html"
    params:
        report_title = "AllTheBacteria's {r} Database Release Creation Report",
        new_db = lambda wildcards: ",".join([f'"{wildcards.o}/allthebacteria-r{wildcards.r}-sigs/allthebacteria-r{wildcards.r}-k{k}.zip"' for k in KSIZES]),
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        scale = config.get('scale_value'),
    conda: "envs/quarto.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell:
        """
        # Mimicing https://github.com/ETH-NEXUS/quarto_example/blob/main/workflow/rules/clean_data_report.smk
        # This will be a stand-alone html document and needs to embed-resources
        # https://quarto.org/docs/output-formats/html-publishing.html#standalone-html
        mkdir -p {wildcards.r}.temp
        cp scripts/report.qmd {wildcards.r}.temp/
        cd {wildcards.r}.temp/

        DIRNAME=$(dirname "{output}")

        quarto render report.qmd \
            -P new_details:{wildcards.r}          -P data_links:{input.data} \
            -P meta_links:{input.meta}            -P new_db:{params.new_db} \
            -P new_mf:{input.new_mf}              -P missing:{input.missing} \
            -P report_title:"{params.report_title}" -P config:"{config}" \
            -P names:{input.names}                -P idents:{input.idents} \
            -P existing:{input.existing}          -P missing_files:{input.missing_files} \
            -P output_dir:{wildcards.o}           -P k_list:{params.k_list}
 
        # the `--output` arg adds an unnecessary `../` to the output file path
        # https://github.com/quarto-dev/quarto-cli/issues/10129
        # just mving it back in place

        mv report.html {output}
        cd ..
        rm -rf {wildcards.r}.temp
        """
