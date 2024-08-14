#! /usr/bin/env python

import sys
import argparse
import pandas as pd

def get_suffix(name):
    ident = name.split(' ')[0]
    assert ident.startswith('GC')
    return ident[3:]

def get_suffix_no_version(name):
    suffix = get_suffix(name)
    assert '.' in suffix
    return suffix.split('.')[0]

def get_float_version(name):
    suffix = get_suffix(name)
    assert '.' in suffix
    return suffix.split('.')[1]

def set_generator(filename):
    with open(filename, 'rt') as fp:
        for i, line in enumerate(fp, start=1):
            if i % 1000 == 0:
                print(f"Reading line {i}", end='\r', flush=True)

            line = line.strip().split('\t')
            if line[0].startswith('#'):
                continue

            accession = line[0]
            yield accession

def update_accessions(metadata_info, update_files):
    """Update accession numbers based on matching suffixes."""

    genbank_set = set(set_generator(update_files[0]))
    genbank_dict = {get_suffix_no_version(item): item for item in genbank_set}
    print('\n', len(genbank_dict))

    def update_row(row):
        suffix = get_suffix_no_version(row['ident'])
        version = get_float_version(row['ident'])
        if suffix in genbank_dict:
            match_item = genbank_dict[suffix]
            new_ver = get_float_version(match_item)
            accession = row['ident'].replace(version, new_ver)
            if version != new_ver:
                print(f"Replacing {row['ident']} with {match_item}")
                print(f'{accession}')
            return accession
        return row['ident']

    metadata_info['ident'] = metadata_info.apply(update_row, axis=1)
    return metadata_info

def main(args):

    metadata_dfs = []
    for metadata_file in args.metadata_files:
        metadata_info = pd.read_csv(metadata_file, header=0, low_memory=False, sep = "\t")
        filtered_metadata_info = metadata_info[["accession", "gtdb_representative", "gtdb_taxonomy"]]
        filtered_metadata_info.loc[:, "accession"] = filtered_metadata_info["accession"].str.replace("RS_", "").str.replace("GB_", "")
        metadata_dfs.append(filtered_metadata_info)

    # Write lineages csv file
    metadata_info = pd.concat(metadata_dfs)
    metadata_info[["superkingdom","phylum","class","order","family","genus","species"]] = metadata_info["gtdb_taxonomy"].str.split(pat=";", expand=True)
    metadata_info.drop(columns=["gtdb_taxonomy"], inplace=True)
    metadata_info.rename(columns={"accession":"ident"}, inplace=True)

    if args.update_files:
        metadata_info = update_accessions(metadata_info, args.update_files)

    metadata_info.to_csv(args.output, sep = ',', index=False)

    if args.reps_csv:
        reps = metadata_info[metadata_info["gtdb_representative"] == 't']
        reps.to_csv(args.reps_csv, sep=',', index=False)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--metadata-files", nargs="+", help="gtdb metadata files")
    p.add_argument("--update-files", nargs="+", help="files generated when updating the database (update_sourmash_dbs.py)")
    p.add_argument("-o", '--output',  help="output lineages csv")
    p.add_argument("-r", '--reps-csv',  help="also output representative lineages csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
