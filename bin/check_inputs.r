#!/usr/bin/env Rscript
# Check the manifest and input samplesheet inputs are complete
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(purrr)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

required <- parser$add_argument_group("Required", "required arguments")

required$add_argument(
  "--input",
  help = "full path to the sample sheet tsv file",
  metavar = "SampleSheet.tsv",
  required = TRUE
)

required$add_argument(
  "--manifest",
  help = "full path to the manifest file",
  metavar = "manifest",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()

if (!file.exists(args$input)) {
  stop("The input samplesheet was not found.")
}

if (!file.exists(args$manifest)) {
  stop("The manifest was not found.")
}

input <- read.delim(args$input)
manifest <- read.delim(args$manifest)

# check manifest paths exist

check_exists <- function(filepath) {
  RCurl::url.exists(filepath) |
  dir.exists(filepath) |
  any(startsWith(filepath, c("gs://", "s3://")))
}

dir_exists <- purrr::pmap_lgl(manifest, ~ check_exists(as.character(..2)))

if (!all(dir_exists)) {
  cat("The following paths were not found: -\n")
  print(manifest[!dir_exists, ])
  stop("Folder paths specified in the manifest were not found.")
} else {
  cat("✓ All paths specified in the manifest were found.\n")
}

# check input samplesheet data present for all keys in manifest
key_in_input <- purrr::map_lgl(
  manifest$key,
  ~ . %in% input$manifest
)

if (!(all(key_in_input))) {
  cat("Input samplesheet data was not found for the following keys: - \n")
  print(manifest[!key_in_input, ]$key)
  stop("Input sample sheet does not contain data for all keys in manifest.")
} else {
  cat("✓ Input samplesheet contains data for all keys in the manifest.\n")
}

cat("Checks passed!\n")

# write the same manifest back out
write.table(manifest,
            "checked_manifest.txt",
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
