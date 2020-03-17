#!/usr/bin/env Rscript
# Finalize SCE with manually revised celltypes
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce_path",
  help = "-path to the SingleCellExperiment",
  metavar = "dir", 
  required = TRUE
)

required$add_argument(
  "--celltype_mappings",
  help = "path to a tsv file with revised celltype mappings",
  metavar = "foo/bar", 
  required = TRUE
)

required$add_argument(
  "--clusters_colname",
  help = "name of the column with cluster numbers",
  metavar = "foo/bar", 
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()

##  ............................................................................
##  Start                                                                   ####

cat(print(tempdir()))

celltype_mappings <- read_celltype_mappings(args$celltype_mappings)

sce <- read_sce(args$sce_path)

sce <- map_custom_celltypes(
    sce, 
    mappings = celltype_mappings, 
    clusters_colname = args$clusters_colname
    )

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "final_sce")
  )

##  ............................................................................
##  Clean up                                                                ####
