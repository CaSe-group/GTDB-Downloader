# GTDB-Downloader
WORK IN PROGRESS Tool to download GTDB samples 
## Description
This tool automates the download of bacterial genomes from the GTDB (Genome Taxonomy Database) via NCBI accession. . It retrieves genome assembly accessions matching a given taxonomic rank and name and gtdb_representative status, then uses the NCBI Datasets CLI tool to download corresponding genomic data. It automatically download (and store) datasets tool and GTDB database. It provides dry run and requre user confirmation before downloading samples.
## Dependendencies 
- R
- data.table R package (should be installed automatically )
- R.utils R package (should be installed automatically )
### Required
- `-n`, `--name_selection`: Taxonomic name to query (e.g., "Escherichia coli")

### Optional
- `-t`, `--taxrank`: Taxonomic rank (default: "species")  
  Options: "phylum", "class", "order", "family", "genus", "species"
- `-r`, `--representative_check`: Filter GTDB representative genomes (default: "t")  
  Options: "t" (yes), "f" (no)
- `-m`, `--domain`: Domain selection (default: "bacteria")  
  Options: "bacteria", "archaea"
- `-db`, `--database`: Custom metadata file path (default: downloads GTDB metadata)
- `-d`, `--dryrun`: Preview without downloading (flag, no value needed)
- `-o`, `--output_dir`: Output directory (default: "output")

## Example Commands
```bash
./fasta_extractor.R -n "Pseudomonas aeruginosa" -t species -r t
```
