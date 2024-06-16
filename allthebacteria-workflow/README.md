A sourmash database for the [AllTheBacteria database](https://doi.org/10.1101/2024.03.08.584059)

The assembly files were all sketched to create a single sourmash database. 

Example release:

- https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/
- https://github.com/AllTheBacteria/AllTheBacteria/releases/tag/v0.2

## Executing 


Create a list of assembly file paths at the ftp server
```
./scripts/ftp_link_list.py ftp.ebi.ac.uk pub/databases/AllTheBacteria/Releases/0.2/assembly -s ftp -p allthebacteria-r0.2-assembly-paths.txt
```

Run the snakemake workflow
```
snakemake -s allthebacteria.smk --use-conda --rerun-incomplete -j 10
```
> A rough estimate of resources required for completing this workflow:
> set to run for 3 days with 10 cpus and 80 GB mem

## Sanity check?

By checking that the total amount of zip files created match the amount expected
```
find allthebacteria-r0.2-sigs/ -maxdepth 2 -type f -name "*.zip" -exec echo "{}" \; | wc -l
```
Returns the count of individual zip files that were created.
For release 0.2, there were `665`.

Which matches `cat allthebacteria-r0.2-assembly-paths.txt | wc -l`


