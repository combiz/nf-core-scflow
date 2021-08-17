#!/usr/bin/env Rscript
# Helper script -- merge multiple tsv into a single tsv
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(purrr)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")

required$add_argument(
  "--filepaths",
  help = "-paths to tsv files",
  metavar = "1.tsv,2.tsv",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$filepaths <- strsplit(args$filepaths, ",")[[1]]

##  ............................................................................
##  Start Merge                                                             ####

tsv_l <- lapply(args$filepaths, read.delim)
merged_tsv <- unique(Reduce(rbind, tsv_l))

write.table(
  merged_tsv,
  file.path(getwd(), "merged.tsv"),
  sep = "\t",
  col.names = TRUE, row.names = FALSE)
