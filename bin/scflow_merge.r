#!/usr/bin/env Rscript
# Merge multiple SingleCellExperiments
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

#options(mc.cores = parallel::detectCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scflow)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce_paths",
  help = "-paths to SingleCellExperiment folders",
  metavar = "dir,dir2", 
  required = TRUE
)

required$add_argument(
  "--ensembl_mappings",
  help = "path to ensembl mappings file",
  metavar = "tsv", 
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$sce_paths <- strsplit(args$sce_paths, ",")[[1]]

##  ............................................................................
##  Start Merge                                                             ####

print(sprintf(
  "Reading %sx SingleCellExperiment's", 
  length(args$sce_paths))
  )

sce_l <- lapply(args$sce_paths, read_sce)

sce <- merge_sce(
  sce_l,
  ensembl_mapping_file = args$ensembl_mappings
  )

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "merged_sce")
  )


##  ............................................................................
##  Clean up                                                                ####

