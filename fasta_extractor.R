#!/usr/bin/env Rscript
require(data.table)

### Parameters
args <- commandArgs(trailingOnly = TRUE)
name_selection <- args[1]
## name_selection <- "Escherichia coli"
output_dir <- "output"
taxrank <- args[2] ## "phylium", "class", "order", "family", "genus", "species"
## taxrank <- "species"
representative_check <- args[3] ## 't' or 'f'
## representative_check <- "t"
dryrun <- args[4] ## "FALSE" "TRUE"
dryrun <- ifelse(length(args) >= 4, as.logical(args[4]), FALSE)

############################
## NOTE: 2. Tool download ##
############################

## File check
file_detection <- list.files()
required_files <- c("bac120_metadata_r220.tsv.gz", "bac120_metadata_r220.tsv", "datasets")
file_check <- required_files %in% file_detection

## GTDB database
if (file_check[1] == FALSE & file_check[2] == FALSE) {
        data <- fread(
                "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz",
                header = TRUE
        )
        fwrite(data, "bac120_metadata_r220.tsv", sep = "\t")
} else if (file_check[1] == TRUE) {
        print("GTDB database found, loading...")
        data <- fread(required_files[1])
} else if (file_check[2] == TRUE) {
        print("GTDB database found, loading...")
        data <- fread(required_files[2])
}

## NCBI dataset
if (file_check[3] == FALSE) {
        download.file("https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets", "datasets")
        system("chmod +x datasets")
} else {
        print("NCBI dataset CLI tool found")
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
        get(taxrank) == name_selection & gtdb_representative == representative_check, c("ncbi_genbank_assembly_accession")
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
fwrite(selection[1:2, ], "download_accession_list.txt", col.names = FALSE)

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
file.remove(rm_list)
system("rm -r ncbi_dataset")

print("Download completed")
