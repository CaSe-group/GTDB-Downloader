# GTDB-Downloader
Tool to download genomes from NCBI based on GTDB taxonomy.

## Description
This tool automates the download of bacterial and archaeal genomes from the GTDB (Genome Taxonomy Database) via NCBI accession. It retrieves genome assembly accessions matching a given taxonomic name, rank, and representative status, then uses the NCBI Datasets CLI tool to download the corresponding genomic data. The NCBI Datasets tool and GTDB database are downloaded and stored automatically. The tool supports dry runs and requires user confirmation before downloading samples.

## Dependencies
- R
- `data.table` R package (installed automatically if missing)
- `R.utils` R package (installed automatically if missing)

## Usage
```bash
./gtdb_downloader.R -n <name> [options]
```

## Parameters

### Required
- `-n`, `--name`: Taxonomic name or identifier to query. Accepts GTDB name, NCBI name, NCBI taxid, or NCBI species taxid (e.g., `"Escherichia coli"`, `"562"`)

### Optional
- `-t`, `--taxrank`: Taxonomic rank to search (default: `"species"`)
  Options: `phylum`, `class`, `order`, `family`, `genus`, `species`
- `-r`, `--representative_check`: Filter for GTDB representative genomes (default: `"t"`)
  Options: `"t"` (yes), `"f"` (no)
- `-m`, `--domain`: Biological domain (default: `"bacteria"`)
  Options: `"bacteria"`, `"archaea"`
- `-o`, `--output_dir`: Output directory (default: `"output"`)
- `--report`: Create a report CSV with accession numbers and taxonomy for each downloaded sample (flag, no value needed)
- `-d`, `--dryrun`: Preview number of matching genomes without downloading (flag, no value needed)
- `-db`, `--database`: Path to a custom local GTDB metadata file
- `-dp`, `--dataset_path`: Path to a pre-installed NCBI Datasets CLI binary (optional; downloaded automatically if not found)
- `-ni`, `--non_interactive`: Skip confirmation prompt before downloading (flag, no value needed)
- `-v`, `--verbose`: Print detailed command output (flag, no value needed)

## Example Commands
```bash
# Download representative genomes for a species
./gtdb_downloader.R -n "Pseudomonas aeruginosa"

# Download all genomes at genus level
./gtdb_downloader.R -n "Pseudomonas" -t genus -r f

# Dry run to preview matches without downloading
./gtdb_downloader.R -n "Escherichia coli" -d

# Query by NCBI taxid with a report
./gtdb_downloader.R -n 562 --report
```

## Acknowledgments
This tool uses metadata from [**GTDB release 232**](https://gtdb.ecogenomic.org/) and the [**NCBI Datasets CLI**](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/) to access the NCBI genome repository. We acknowledge the valuable contributions of the GTDB and NCBI teams in providing these resources.
