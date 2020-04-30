#!/usr/bin/env Rscript
#' Annotate integrated, reduced dimension,
#' and clustered SingleCellExperiment object
# Mahdi Moradi Marjaneh

# ____________________________________________________________________________
# Initialization ####

options(mc.cores = future::availableCores())

## ............................................................................
## Load packages ####
library(argparse)
library(scFlow)
library(parallel)

## ............................................................................
## Parse command-line arguments ####

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
  "--categorical_covariates",
  help = "-a list of categorical variables",
  metavar = "individual,diagnosis,region,sex",
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args ####

args <- parser$parse_args()
args <- purrr::map(args, function(x) {
if (length(x) == 1) {
if (toupper(x) == "TRUE") {
return(TRUE)
}
if (toupper(x) == "FALSE") {
return(FALSE)
}
if (toupper(x) == "NULL") {
return(NULL)
}
}
return(x)
})

## ............................................................................
## Annotate integrated sce ####

sce <- read_sce(args$sce_path)

sce <- annotate_integrated_sce(
sce,
categorical_covariates = args$categorical_covariates
)

## ............................................................................
## Save Outputs ####

# Save SingleCellExperiment
write_sce(
sce = sce,
folder_path = file.path(getwd(), "integrated_sce")
)