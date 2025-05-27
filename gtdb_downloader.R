#!/usr/bin/env Rscript
suppressMessages(
        for (package in c("data.table", "R.utils")) {
                if (!require(package, character.only = T, quietly = T)) {
                        install.packages(package)
                        library(package, character.only = T)
                }
        }
)

#################
## Update info ##
#################
## In case of database update please modify just the line annotated with following text
## (if no column or parameters were modified):
# NOTE: DATABASE UPDATE MODIFY JUST THIS LINES

### Parameters
args <- commandArgs(trailingOnly = TRUE)

### Defaults
name_selection <- "623"
taxrank <- "species"
representative <- "t"
dryrun <- FALSE
output_dir <- "output"
domain <- "bacteria"
database_file <- NULL
non_interactive <- FALSE
report <- FALSE
verbose <- FALSE

##################
## NOTE: 1. SRC ##
##################
green <- function(x) paste0("\033[1;32m", x, "\033[0m")
yellow <- function(x) paste0("\033[1;33m", x, "\033[0m")
red <- function(x) paste0("\033[1;31m", x, "\033[0m")

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

tax_split_ncbi <- function(dt) {
        ## tax split
        tax_levels <- c("ncbi_phylium", "ncbi_class", "ncbi_order", "ncbi_family", "ncbi_genus", "ncbi_species")
        dt[, (tax_levels) := tstrsplit(ncbi_taxonomy, ";")[2:7]]
        ## Pref remove
        for (col in tax_levels) {
                dt[, (col) := sub("^[a-z]__", "", get(col))]
        }
        return(dt)
}


# Parse parameters arguments
i <- 1
while (i <= length(args)) {
        if (args[i] %in% c("--name", "--name_selection", "-n")) {
                name_selection <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--taxrank", "-t")) {
                taxrank <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--representative_check", "-r")) {
                representative <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--dryrun", "-d")) {
                dryrun <- TRUE
                i <- i + 1
        } else if (args[i] %in% c("--output_dir", "-o")) {
                output_dir <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--domain", "-m")) {
                domain <- tolower(args[i + 1])
                if (!domain %in% c("bacteria", "archaea")) {
                        stop("Domain must be either 'bacteria' or 'archaea'")
                }
                i <- i + 2
        } else if (args[i] %in% c("--database", "-db")) {
                database_file <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--non_interactive", "-ni")) {
                non_interactive <- TRUE
                i <- i + 1
        } else if (args[i] %in% c("--report")) {
                report <- TRUE
                i <- i + 1
        } else if (args[i] %in% c("--dataset_path", "-dp")) {
                dataset_path <- args[i + 1]
                i <- i + 2
        } else if (args[i] %in% c("--verbose", "-v")) {
                verbose <- TRUE
                i <- i + 1
        } else if (args[i] %in% c("--help", "-h")) {
                cat("Usage: fasta_extractor [options]\n")
                cat("Example command:\n")
                cat('./gtdb_downloader.R -n "Pseudomonas aeruginosa" \n')
                cat("Options:\n")
                cat("  -n, --name <name>   Name of the taxon to select based on GTDB (preferable), NCBI taxonomy or NCBI taxid (required; eg. Escherichia coli)\n")
                cat("  -t, --taxrank <rank>         Taxonomic rank [default: species]\n")
                cat("  -r, --representative_check <t/f> Check for representative genomes [default: t]\n")
                cat("  -d, --dryrun                 Dry run (no actual processing) [default: FALSE]\n")
                cat("  -o, --output_dir <dir>       Output directory [default: output]\n")
                cat("  -m, --domain <domain>        Domain (bacteria or archaea) [default: bacteria]\n")
                cat("  -db, --database <file>       Custom database file path\n")
                cat("  -dp, --dataset_path          Path to dataset tool (not required)\n")
                cat("  -ni, --non_interactive       Do not ask for confirmation before downloading genomes\n")
                cat("  --report                     Create a report with each sample accession number and taxonomy\n")
                cat("  -v, --verbose                Verbose command output\n")
                cat("  -h, --help                   Show this help message\n")
                quit()
        } else {
                # Fallback to positional arguments
                if (is.null(name_selection)) {
                        name_selection <- args[i]
                } else if (i == 2) {
                        taxrank <- args[i]
                } else if (i == 3) {
                        representative <- args[i]
                } else if (i == 4) dryrun <- as.logical(args[i])
                i <- i + 1
        }
}

# Check if name_selection was given
if (is.null(name_selection)) {
        stop(red("--name argument is required\nUse --help for usage information"), call. = FALSE)
}
### Parameters
cat("Parameters:", "\n")
cat("Name:", name_selection, "\n")
cat("Representative mode:", representative, "\n")
cat("Taxonomical rank:", taxrank, "\n")
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
                default_files <- c("bac120_metadata_r226.tsv.gz", "bac120_metadata_r226.tsv")
                download_url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_metadata_r226.tsv.gz"
        } else if (domain == "archaea") {
                default_files <- c("ar53_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv")
                download_url <- "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_metadata_r226.tsv.gz"
        }
        ## File check
        file_detection <- list.files()
        file_check <- default_files %in% file_detection

        if (!any(file_check)) {
                cat("Downloading", domain, "metadata from GTDB...\n")
                ## Trimming the database
                data <- fread(download_url, header = TRUE)
                cat("Reducing the database footprint...\n")
                col_to_preserve <- c(
                        "gtdb_taxonomy", "ncbi_taxonomy",
                        "gtdb_representative", "ncbi_genbank_assembly_accession",
                        "ncbi_taxid", "ncbi_species_taxid"
                )
                data <- subset(data, , col_to_preserve)
                data <- tax_split(data)
                data <- tax_split_ncbi(data)
                local_file <- default_files[2] # Save as uncompressed TSV
                fwrite(data, local_file, sep = "\t")
                cat("Saved database to:", local_file, "\n")
        } else {
                existing_file <- default_files[which(file_check)[1]]
                cat(green(paste("GTDB", domain, "database found, loading:", existing_file, "\n")))
                data <- fread(existing_file)
        }
}

## NCBI dataset
if (!file.exists("datasets") & !exists(quote(dataset_path))) {
        cat("Downloading NCBI datasets CLI tool...\n")
        download.file(
                "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets",
                "datasets",
                quiet = TRUE
        )
        system("chmod +x datasets", ignore.stdout = TRUE, ignore.stderr = TRUE)
        cat(green("NCBI datasets CLI tool installed and made executable\n"))
} else {
        cat(green("NCBI datasets CLI tool found\n"))
}

########################
## NOTE: 3. File input##
########################
## List prep
if (!("species" %in% colnames(data)) && !("ncbi_species" %in% colnames(data))) {
        processed_data <- tax_split(data)
        processed_data <- tax_split_ncbi(processed_data)
} else {
        processed_data <- data
        rm(data)
}
## Tax id
## Selection fallback taxrank > ncbi_taxid > ncbi_species_taxid

selection <- processed_data[
        get(taxrank) == name_selection & gtdb_representative == representative,
        c("ncbi_genbank_assembly_accession")
]

if (nrow(selection) == 0) {
        selection <- processed_data[
                paste0("ncbi_", get(taxrank)) == name_selection & gtdb_representative == representative,
                c("ncbi_genbank_assembly_accession")
        ]
        if (verbose == TRUE) {
                cat(yellow("No GTDB taxonomy name found falling back to NCBI taxonomical rank name\n"))
                cat(yellow("For the optimal results please use GTDB taxonomy ranks\n"))
        }
}

if (nrow(selection) == 0) {
        if (verbose == TRUE) {
                cat(yellow("No NCBI Species Taxid found falling back to NCBI Taxid\n"))
        }
        if (representative == "f") {
                ncbi_name <- processed_data[ncbi_taxid == name_selection, ncbi_species][1]
                selection <- processed_data[
                        ncbi_taxid == name_selection,
                        c("ncbi_genbank_assembly_accession")
                ]
        } else if (representative == "t") {
                ncbi_name <- processed_data[ncbi_taxid == name_selection, ncbi_species][1]
                selection <- processed_data[
                        species == ncbi_name &
                                gtdb_representative == representative,
                        c("ncbi_genbank_assembly_accession")
                ]
        }
        if (verbose == TRUE) {
                cat(yellow("Taxid used with GTDB representative genome...\n"))
        }

        if (nrow(selection) > 0) {
                cat(yellow(paste0(
                        "Guessing GTDB species using NCBI taxonomy\n",
                        "Guessed species: ", ncbi_name, "\n"
                )))
        }
}

if (nrow(selection) == 0) {
        if (verbose == TRUE) {
                cat(yellow("No NCBI Taxid found falling back to Species Taxid \n"))
        }
        if (representative == "f") {
                selection <- processed_data[
                        ncbi_species_taxid == name_selection,
                        c("ncbi_genbank_assembly_accession")
                ]
        } else if (representative == "t") {
                ncbi_name <- processed_data[ncbi_species_taxid == name_selection, ncbi_species][1]
                selection <- processed_data[
                        species == ncbi_name & gtdb_representative == representative,
                        c("ncbi_genbank_assembly_accession")
                ]
                if (verbose == TRUE) {
                        cat(yellow("Taxid used with GTDB representative genome...\n"))
                }
                if (nrow(selection) > 0) {
                        cat(yellow(paste0(
                                "Guessing GTDB species using NCBI taxonomy\n",
                                "Guessed species:", ncbi_name, "\n"
                        )))
                }
        }
}

if (nrow(selection) == 0 & representative == "t") {
        ncbi_name <- processed_data[ncbi_species_taxid == name_selection, ncbi_species][1]

        gtdb_name <- processed_data[
                ncbi_species == ncbi_name,
                c("species")
        ]
        species_count <- as.data.frame(table(gtdb_name))
        if (nrow(species_count) > 1 & nrow(species_count != 0)) {
                setcolorder(processed_data, c("gtdb_taxonomy", "ncbi_taxonomy"), after = ncol(processed_data))
                setcolorder(processed_data, "ncbi_genbank_assembly_accession", before = 1)
                fwrite(processed_data[processed_data$ncbi_genbank_assembly_accession %in% selection$ncbi_genbank_assembly_accession, ], file = paste0(output_dir, "/ambiguise_samples_report.csv"))
                cat(red(paste0("Ambiguise samples found based on pure NCBI taxonomy\n", "Sample count found: ", nrow(selection), "\n Ambiguise samples output printed to:", output_dir, "/ambiguise_samples_report.csv file")))
                quit()
        } else {
                selection <- processed_data[
                        species == gtdb_name[1] & gtdb_representative == "t",
                        c("ncbi_genbank_assembly_accession")
                ]
        }

        cat(green(paste0(
                "Using pure NCBI taxonomy: ", "\n",
                "NCBI Species:", ncbi_name, "\n",
                "GTDB Species:", gtdb_name[1], "\n"
        )))
}

## User confirmation and sample count check
n_records <- nrow(selection)

if (n_records == 0) {
        ## Print similarity matrix if n_records are 0
        cat(red(paste0("\n", name_selection, " ", n_records, " records found.\n")))
        unique_species <- processed_data[, .N, by = get(taxrank)][order(-N)]
        colnames(unique_species) <- c("Name", "N")
        unique_species[, distance := adist(Name, name_selection)]
        setorder(unique_species, distance)
        top_similar <- unique_species[Name != name_selection][1:min(10, .N)]
        has_similar <- nrow(top_similar) > 0
        if (has_similar) {
                cat("\nSimilar records:\n(Not selected for download)\n")
                print(top_similar[, .(Name, Records = N)])
        } else {
                cat(red("No similar species found\n"))
        }
        quit(status = 0)
} else {
        cat(green(paste0("\n", name_selection, " ", n_records, " records found.\n")))
}
if (dryrun == TRUE) {
        cat(yellow("Dryrun mode activated quitting"))
        quit()
} else {
        cat("Would you like to proceed with the download? [y/n]: ")
}
if (non_interactive == TRUE) {
        print("Non interactive mode say yes")
} else {
        response <- readLines("stdin", n = 1)
        if (!tolower(response) %in% c("y", "yes")) {
                cat("Aborting at user request.\n")
                quit(status = 0)
        }
}

## Accession write (for datasets tool)
fwrite(selection, "download_accession_list.txt", col.names = FALSE)

## Genomes download
if (!exists(quote(dataset_path))) {
        dataset_path_input <- "./datasets "
} else {
        dataset_path_input <- paste0(dataset_path, " ")
}
if (verbose == TRUE) {
        invisible(system(paste0(dataset_path_input, "download genome accession --inputfile download_accession_list.txt --dehydrated --include genome")))
} else {
        system(paste0(dataset_path_input, "download genome accession --inputfile download_accession_list.txt --dehydrated --include genome"))
}


unzip("ncbi_dataset.zip")

if (verbose == TRUE) {
        system(paste0(dataset_path_input, "rehydrate --directory ."))
} else {
        invisible(system(paste0(dataset_path_input, "rehydrate --directory .")))
}

dir.create(output_dir)
invisible(file.rename(
        from = list.files(path = "ncbi_dataset/data/", pattern = "\\.fna$", recursive = TRUE, full.names = TRUE),
        to = file.path(paste0("./", output_dir), basename(list.files("ncbi_dataset/data/", "\\.fna$", recursive = TRUE)))
))

if (report == TRUE) {
        setcolorder(processed_data, c("gtdb_taxonomy", "ncbi_taxonomy"), after = ncol(processed_data))
        setcolorder(processed_data, "ncbi_genbank_assembly_accession", before = 1)

        fwrite(processed_data[processed_data$ncbi_genbank_assembly_accession %in% selection$ncbi_genbank_assembly_accession, ], file = paste0(output_dir, "/download_report.csv"))
}

## Cleanup
rm_list <- c("download_accession_list.txt", "ncbi_dataset.zip", "md5sum.txt", "README.md")
invisible(file.remove(rm_list))
invisible(system("rm -r ncbi_dataset"))

cat(green("Download completed"))
