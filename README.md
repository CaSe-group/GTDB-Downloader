# GTDB-Downloader
WORK IN PROGRESS Tool to download GTDB samples 
## Description
This tool automates the download of bacterial genomes from the GTDB (Genome Taxonomy Database) via NCBI accession. . It retrieves genome assembly accessions matching a given taxonomic rank and name and gtdb_representative status, then uses the NCBI Datasets CLI tool to download corresponding genomic data. It automatically download (and store) datasets tool and GTDB database. It provides dry run and requre user confirmation before downloading samples.
## Dependendencies 
- R
- data.table R package (should be installed automatically (I hope))
## Input
- Name: Taxonomic name to query (e.g., "Escherichia coli").
- Taxonomic rank: Taxonomic rank to filter by (options: "phylium", "class", "order", "family", "genus", "species").
- Representative: Filter GTDB representative genomes ("t" for yes, "f" for no).
- Dryrun (optional): Preview results without downloading - dryrun FALSE still require confirmation before downloading ("TRUE" or "FALSE").
## Example command
``` bash
./fasta_extractor.R 'Pseudomonas aeruginosa' 'species' 't'
``` 
