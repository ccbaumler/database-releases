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
            selected_line = random.choice(file_list)
            dir_name_list = [os.path.basename(selected_line).split(".")[0]]

            directory_path = f"test/allthebacteria-r{release}"

            if not os.path.exists(directory_path):
                os.makedirs(directory_path)
            else:
                print(f"The directory '{directory_path}' already exists.")
    
            test_filename = f"{directory_path}/sub-{file}"
            
            with open(test_filename, 'w') as fp:
                fp.write(f"{selected_line}")

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

rule check_sketch:
    input:
        zip_file = "{o}/allthebacteria-r{r}-sigs/{dir_name}/{dir_name}.zip",
    log:
        '{o}/logs/check_sketch/allthebacteria-r{r}-sig-{dir_name}.log'
    shell:"""
        unzip -v {input.zip_file} 2>&1 | tee -a {log}
    """

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
    priority: 1
    shell: """
        {input.script} ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/{wildcards.r}/metadata -s ftp -p {output.meta}

        {input.script} ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/{wildcards.r}/assembly -s ftp -p {output.data}
    """

# create a 1 sub_assembly file for testing!!!
#rule test_with_sub_assembly_summary:
#    input:
#        data = "{o}/allthebacteria-r{r}-data-links.txt",
#    output:
#        data = "{o}/allthebacteria-r{r}-sub-data-links.txt",
#    shell:"""
#        echo {input.data}
#        cat {input.data} | wc -l
#
##        awk 'BEGIN {{srand();}} {{print rand() " " $0}}' "{input.data}" | sort -n | cut -d" " -f 2 | tail --lines=1 > "{output.data}"
#
#        cat {output.data} | wc -l
#        cat {output.data}
#    """

rule download_links:
    input:
        unpack(getInputFilesForLinks),
        meta = "{o}/allthebacteria-r{r}-metadata-links.txt",
    output:
        tsv = "{o}/allthebacteria-r{r}-metadata/sylph.tsv.gz",
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
	partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    priority: 1
    threads: 
        10
    shell: """
        # read the txt files, use only 1 link (-n), parallelize wget cmd upto 10 times (-P), path for downloads (-P), if interrupted continue the download (--continue), do not download anything that exists (--no-clobber), try up to 3 times per download (--tries), wait 1 sec per try (--wait), quiet output (-nv) 
        cat {input.data} | xargs -n 1 -P {threads} wget -P {wildcards.o}/allthebacteria-r{wildcards.r}-data/ --continue --no-clobber --tries=3 --wait=1 -nv && cat {input.meta} | xargs -n 1 -P {threads} wget -P {wildcards.o}/allthebacteria-r{wildcards.r}-metadata/ --continue --no-clobber --tries=3 --wait=1 -nv
    """

rule extract_info:
    input:  
        info = "{o}/allthebacteria-r{r}-metadata/sylph.tsv.gz",
    output:
        info = "{o}/allthebacteria-r{r}-metadata/sylph.tsv",
    resources:
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    priority: 1
    shell:
        """
        gzip -kvdf {input.info}
    """

rule extract_files:
    input:
        placeholder = "{o}/allthebacteria-r{r}-metadata/sylph.tsv",
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
    shell:
        """
        echo "Extracting sequence files from compressed archives"
        pv {params.data} | tar --use-compress-program="xz -T0 -q" --skip-old-files -xf - -C {params.outdir}
        #tar --skip-old-files -xvf {params.data} -C {params.outdir}
        echo "Extracted $(ls -1 {output.seq_dir} | wc -l) files!"
        """

#info_dict = None
#
#def load_info_dict(info_files=[]):
#    global info_dict
#    if info_dict is None:
#        info_dict = {}
#        for info_file in info_files:
#            with open(info_file, 'r') as fp:
#                reader = csv.DictReader(fp, delimiter='\t')
#                for row in reader:
#                    info_dict[row['sample']] = row['Contig_name'] + " " + row['sample']
#
#load_info_dict([f"allthebacteria-r{release}-metadata/sylph.tsv" for release in RELEASES],)

rule build_csv:
    input:
        seq_dir = "{o}/allthebacteria-r{r}-seqs/{dir_name}",
        info_file = "{o}/allthebacteria-r{r}-metadata/sylph.tsv"
    output:
        csv_file = "{o}/allthebacteria-r{r}-sigs/{dir_name}/manysketch.csv"
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    #group: "group"
    params:
        dir_name = "{dir_name}"
    priority: 1
    run:

        def load_info_dict(info_files=[]):
            info_dict = {}
            for info_file in info_files:
                with open(info_file, 'r') as fp:
                    reader = csv.DictReader(fp, delimiter='\t')
                    for row in reader:
                        info_dict[row['sample']] = row['Contig_name'] + " " + row['sample']
            return info_dict
        #info_dict = {}

        info_dict = load_info_dict([f"{wildcards.o}/allthebacteria-r{release}-metadata/sylph.tsv" for release in RELEASES],)

        with open(output.csv_file, 'w', newline='') as csvfile:
             csvwriter = csv.writer(csvfile)
             csvwriter.writerow(['name', 'genome_filename', 'protein_filename'])

             filepaths = [os.path.join(input.seq_dir, f) for f in os.listdir(input.seq_dir) if os.path.isfile(os.path.join(input.seq_dir, f))]
             num_files = len(filepaths)

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
        allowed_jobs=10, #lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    threads: 16
    conda: 'envs/branchwater.yaml'
    priority: 1
    params:
        k_list = lambda wildcards: ",".join([f"k={k}" for k in KSIZES]),
        scale = {SCALE},
    shell:"""
        sourmash scripts manysketch {input.csv_file} -p {params.k_list},scaled={params.scale},abund -o {output.zip_file}
    """

rule cat_all:
    input:
        zip_file = expand("{o}/allthebacteria-r{r}-sigs/{dir_name}/{dir_name}.zip", o=OUTDIR, r=RELEASES, dir_name=[dn for r in RELEASES for dn in dir_name_dict[r]]),
    output:
        all = "{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}.zip",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    conda: 'envs/sourmash.yaml',
    params:
        path = "{o}/allthebacteria-r{r}-sigs"
    shell:"""
        find {params.path} \
        -maxdepth 2 -type f -name "*.zip" -exec sh -c \
        'sourmash sig cat "$@" -o "$0" ' "{output.all}" {{}} +
    """

rule cat_by_k:
    input:
        zip_file = "{o}/allthebacteria-r{r}-sigs/allthebacteria-r{r}.zip",
    output:
        k_zip_file = "{o}/allthebacteria-r{r}-k{k}.zip",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs = 100,
        partition = "bmh",
    conda: 'envs/sourmash.yaml',
    shell:"""
        sourmash sig cat {input.zip_file} -k {wildcards.k} -o {output.k_zip_file}
    """
