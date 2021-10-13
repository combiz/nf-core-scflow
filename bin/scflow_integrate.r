#!/usr/bin/env Rscript
# Integrate multiple single cell datasets (samples)
# Mahdi Moradi Marjaneh

# ____________________________________________________________________________
# Initialization ####

options(mc.cores = max(2, future::availableCores(methods = "mc.cores")))

## ............................................................................
## Load packages ####
library(scFlow)
library(argparse)
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
  "--method",
  required = TRUE,
  help ="The integration method to use",
  metavar = "Liger"
)

required$add_argument(
  "--unique_id_var",
  required = TRUE,
  help ="Unique id variable",
  metavar = "manifest"
)

required$add_argument(
  "--take_gene_union",
  default = FALSE,
  required = TRUE,
  help ="Whether to fill out raw.data matrices with union of genes across all datasets (filling in 0 for missing data)",
  metavar = "Boolean"
)

required$add_argument(
  "--remove_missing",
  default = TRUE,
  required = TRUE,
  help ="Whether to remove cells not expressing any measured genes, and genes not expressed in any cells",
  metavar = "Boolean"
)

required$add_argument(
  "--num_genes",
  default = 3000,
  type = "integer",
  required = TRUE,
  help ="Number of genes to find for each dataset",
  metavar = "N"
)

required$add_argument(
  "--combine",
  default = "union",
  required = TRUE,
  help ="How to combine variable genes across experiments",
  metavar = "union,intersect"
)

required$add_argument(
  "--capitalize",
  default = FALSE,
  required = TRUE,
  help ="Capitalize gene names to match homologous genes(ie. across species)",
  metavar = "Boolean"
)

required$add_argument(
  "--use_cols",
  default = TRUE,
  required = TRUE,
  help ="Treat each column as a cell",
  metavar = "Boolean"
)

required$add_argument(
  "--k",
  default = 30,
  type = "integer",
  required = TRUE,
  help ="Inner dimension of factorization (number of factors)",
  metavar = "N"
)

required$add_argument(
  "--lambda",
  default = 5.0,
  type = "double",
  required = TRUE,
  help ="Regularization parameter. Larger values penalize dataset-specific effects more strongly (ie. alignment should increase as lambda increases)",
  metavar = "N"
)

required$add_argument(
  "--thresh",
  default = 0.0001,
  type = "double",
  required = TRUE,
  help ="Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh",
  metavar = "N"
)

required$add_argument(
  "--max_iters",
  default = 100,
  type = "integer",
  required = TRUE,
  help ="Maximum number of block coordinate descent iterations to perform",
  metavar = "N"
)

required$add_argument(
  "--nrep",
  default = 1,
  type = "integer",
  required = TRUE,
  help ="Number of restarts to perform",
  metavar = "N"
)

required$add_argument(
  "--rand_seed",
  default = 1,
  type = "integer",
  required = TRUE,
  help ="Random seed to allow reproducible results",
  metavar = "N"
)

required$add_argument(
  "--knn_k",
  default = 20,
  type = "integer",
  required = TRUE,
  help ="Number of nearest neighbors for within-dataset knn graph",
  metavar = "N"
)

required$add_argument(
  "--ref_dataset",
  default = '',
  required = TRUE,
  help ="Name of dataset to use as a reference for normalization",
  metavar = "ref"
)

required$add_argument(
  "--min_cells",
  default = 2,
  type = "integer",
  required = TRUE,
  help ="Minimum number of cells to consider a cluster shared across datasets",
  metavar = "N"
)

required$add_argument(
  "--quantiles",
  default = 50,
  type = "integer",
  required = TRUE,
  help ="Number of quantiles to use for quantile normalization",
  metavar = "N"
)

required$add_argument(
  "--resolution",
  default = 1,
  type = "double",
  required = TRUE,
  help ="Controls the number of communities detected (Higher resolution -> more communities)",
  metavar = "N"
)

required$add_argument(
  "--center",
  default = FALSE,
  required = TRUE,
  help ="Centers the data when scaling factors (useful for less sparse modalities like methylation data)",
  metavar = "Boolean"
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
## Integrate sce ####

sce <- read_sce(args$sce_path)

sce <- integrate_sce(
  sce,
  method = args$method,
  unique_id_var = args$unique_id_var,
  take_gene_union = args$take_gene_union,
  remove.missing = args$remove_missing,
  num_genes = args$num_genes,
  combine = args$combine,
  capitalize = args$capitalize,
  use_cols = args$use_cols,
  num_cores = future::availableCores(methods = "mc.cores"),
  k = args$k,
  lambda = args$lambda,
  thresh = args$thresh,
  max_iters = args$max_iters,
  nrep = args$nrep,
  H_init = NULL,
  W_init = NULL,
  V_init = NULL,
  rand_seed = args$rand_seed,
  knn_k = args$knn_k,
  ref_dataset = args$ref_dataset,
  min_cells = args$min_cells,
  quantiles = args$quantiles,
  resolution = args$resolution,
  center = args$center,
  print_obj = FALSE
)


## ............................................................................
## Save Outputs ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "integrated_sce"),
  write_metadata = TRUE
)
