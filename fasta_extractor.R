#!/usr/bin/env Rscript

suppressMessages(require(data.table))
suppressMessages(require(R.utils))

### Parameters
args <- commandArgs(trailingOnly = TRUE)

### Defaults 
name_selection <- NULL
taxrank <- "species"
representative <- "t"
dryrun <- FALSE
output_dir <- "output"
domain <- "bacteria"  
database_file <- NULL 

# Parse named arguments
i <- 1
while (i <= length(args)) {
  if (args[i] %in% c("--name_selection", "-n")) {
    name_selection <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--taxrank", "-t")) {
    taxrank <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--representative_check", "-r")) {
    representative_check <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--dryrun", "-d")) {
    dryrun <- TRUE
    i <- i + 1
  } else if (args[i] %in% c("--output_dir", "-o")) {
    output_dir <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--domain", "-m")) {  
    domain <- tolower(args[i+1])
    if (!domain %in% c("bacteria", "archaea")) {
      stop("Domain must be either 'bacteria' or 'archaea'")
    }
    i <- i + 2
  } else if (args[i] %in% c("--database", "-db")) {  
    database_file <- args[i+1]
    i <- i + 2
  } else if (args[i] %in% c("--help", "-h")) {
    cat("Usage: fasta_extractor [options]\n")
    cat("Options:\n")
    cat("  -n, --name_selection <name>   Name of the taxon to select (required)\n")
    cat("  -t, --taxrank <rank>         Taxonomic rank [default: species]\n")
    cat("  -r, --representative_check <t/f> Check for representative genomes [default: t]\n")
    cat("  -d, --dryrun                 Dry run (no actual processing) [default: FALSE]\n")
    cat("  -o, --output_dir <dir>       Output directory [default: output]\n")
    cat("  -m, --domain <domain>        Domain (bacteria or archaea) [default: bacteria]\n")
    cat("  -db, --database <file>       Custom database file path\n")
    cat("  -h, --help                   Show this help message\n")
    quit()
  } else {
    # Fallback to positional arguments
    if (is.null(name_selection)) name_selection <- args[i]
    else if (i == 2) taxrank <- args[i]
    else if (i == 3) representative_check <- args[i]
    else if (i == 4) dryrun <- as.logical(args[i])
    i <- i + 1
  }
}

# Check if name_selection was givven
if (is.null(name_selection)) {
  stop("Error: --name_selection argument is required\nUse --help for usage information", call. = FALSE)
}
############################
## NOTE: 2. Tool download ##
############################

### Database
if (!is.null(database_file)) {
  # Custom database
  if (!file.exists(database_file)) {
    stop(paste("Specified database file not found:", database_file))
  }
  cat("Loading custom database from:", database_file, "\n")
  data <- fread(database_file)
} else {
  # NOTE: DATABASE UPDATE MODIFY JUST THIS LINES
  # Default GTDB database
  if (domain == "bacteria") {
    default_files <- c("bac120_metadata_r220.tsv.gz", "bac120_metadata_r220.tsv")
    download_url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
  } else if (domain == "archaea") {
    default_files <- c("ar53_metadata_r220.tsv.gz", "ar53_metadata_r220.tsv")
    download_url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/ar53_metadata_r220.tsv.gz"
  }

## File check
  file_detection <- list.files()
  file_check <- default_files %in% file_detection

  if (!any(file_check)) {
    cat("Downloading", domain, "metadata from GTDB...\n")
    data <- fread(download_url, header = TRUE)
    local_file <- default_files[2]  # Save as uncompressed TSV
    fwrite(data, local_file, sep = "\t")
    cat("Saved database to:", local_file, "\n")
  } else {
    existing_file <- default_files[which(file_check)[1]]
    cat("GTDB", domain, "database found, loading:", existing_file, "\n")
    data <- fread(existing_file)
  }
}

## NCBI dataset
if (!file.exists("datasets")) {
  cat("Downloading NCBI datasets CLI tool...\n")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets",
    "datasets",
    quiet = TRUE  
  )
  system("chmod +x datasets", ignore.stdout = TRUE, ignore.stderr = TRUE)
  cat("NCBI datasets CLI tool installed and made executable\n")
} else {
  cat("NCBI datasets CLI tool found\n")
}

#################
## NOTE: 2. SRC##
#################

tax_split <- function(dt) {
        ## tax split
        tax_levels <- c("phylium", "class", "order", "family", "genus", "species")
        dt[, (tax_levels) := tstrsplit(gtdb_taxonomy, ";")[2:7]]
        ## Pref remove
        for (col in tax_levels) {
                dt[, (col) := sub("^[a-z]__", "", get(col))]
        }
        return(dt)
}

########################
## NOTE: 3. File input##
########################
## List prep
processed_data <- tax_split(data)

selection <- processed_data[
        get(taxrank) == name_selection & gtdb_representative == representative, c("ncbi_genbank_assembly_accession")
]

## User confirmation and sample count check
n_records <- nrow(selection)
unique_species <- processed_data[, .N, by = get(taxrank)][order(-N)]
colnames(unique_species) <- c("Name", "N")
unique_species[, distance := adist(Name, name_selection)]
setorder(unique_species, distance)
top_similar <- unique_species[Name != name_selection][1:min(10, .N)]

if (dryrun == TRUE) {
        cat(paste0("\n", name_selection, " ", n_records, " records found.\n"))
        if (nrow(top_similar) > 0) {
                cat("\nSimilar records:\n(Not selected for download)\n")
                print(top_similar[, .(Name = Name, Records = N)])
                quit(status = 0)
        } else {
                cat("\nNo similar species found\n")
        }
}

cat(paste0("\n", name_selection, " ", n_records, " records found.\n"))
if (nrow(top_similar) > 0) {
        cat("\nSimilar records:\n(Not selected for download)\n")
        print(top_similar[, .(Name = Name, Records = N)])
} else {
        cat("\nNo similar species found\n")
}
cat("Would you like to proceed with the download? [y/n]: ")
response <- readLines("stdin", n = 1)
if (!tolower(response) %in% c("y", "yes")) {
        cat("Aborting at user request.\n")
        quit(status = 0)
}

## Accession write (for datasets tool)
fwrite(selection, "download_accession_list.txt", col.names = FALSE)

## Genomes download
system("./datasets download genome accession --inputfile download_accession_list.txt --dehydrated --include genome")
unzip("ncbi_dataset.zip")
system("./datasets rehydrate --directory .")
dir.create(output_dir)
file.rename(
        from = list.files(path = "ncbi_dataset/data/", pattern = "\\.fna$", recursive = TRUE, full.names = TRUE),
        to = file.path(paste0("./", output_dir), basename(list.files("ncbi_dataset/data/", "\\.fna$", recursive = TRUE)))
)

## Cleanup
rm_list <- c("download_accession_list.txt", "ncbi_dataset.zip", "md5sum.txt", "README.md")
suppressMessages(file.remove(rm_list))
suppressMessages(system("rm -r ncbi_dataset"))

print("Download completed")
