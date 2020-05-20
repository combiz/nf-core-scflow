#!/usr/bin/env Rscript
# Reduce dimensions for a SCE
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(parallel)
library(SingleCellExperiment) # due to monocle3 missing namespace::

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
  "--input_reduced_dim",
  help = "input reducedDim to use for further dim reds",
  metavar = "PCA,Liger", 
  required = TRUE
)

required$add_argument(
  "--reduction_methods",
  help = "methods to use for dimensionality reduction",
  metavar = "PCA,tSNE,UMAP,UMAP3D", 
  required = TRUE
)

required$add_argument(
  "--vars_to_regress_out",
  help = "variables to regress out before finding singlets",
  metavar = "nCount_RNA,pc_mito", 
  required = TRUE
)

required$add_argument(
  "--pca_dims",
  type = "integer", 
  default = 20,
  help = "the number of PCA dimensions used",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--n_neighbors",
  type = "integer", 
  default = 30,
  help = "the number of nearest neighbors",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--n_components",
  type = "integer", 
  default = 2,
  help = "the number of UMAP dimensions (2 or 3)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--init",
  help = "type of initialization for UMAP coordinates (uwot)",
  metavar = "pca", 
  required = TRUE
)

required$add_argument(
  "--metric",
  help = "type of distance metric for nearest neighbours (uwot)",
  metavar = "euclidean", 
  required = TRUE
)

required$add_argument(
  "--n_epochs",
  type = "integer", 
  default = 500,
  help = "number of epochs for optimization of embeddings (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--learning_rate",
  type = "double", 
  default = 1.0,
  help = "initial learning rate used in optimization (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--min_dist",
  type = "double", 
  default = 0.3,
  help = "effective minimum distance between embedded points (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--spread",
  type = "double", 
  default = 1.0,
  help = "effective scale of embedded points (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--set_op_mix_ratio",
  type = "double", 
  default = 1.0,
  help = "interpolation between fuzzy union and intersection set operation (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--local_connectivity",
  type = "integer", 
  default = 1,
  help = "number of nearest neighbours assumed connected locally (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--repulsion_strength",
  type = "double", 
  default = 1.0,
  help = "weighting applied to negative samples in optimization (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--negative_sample_rate",
  type = "double", 
  default = 5.0,
  help = "number of negative edge samples per positive edge (uwot)",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--fast_sgd",
  help = "faster but less reproducible UMAP (uwot) (lgl)",
  metavar = "FALSE", 
  required = TRUE
)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$fast_sgd <- as.logical(args$fast_sgd)
args$input_reduced_dim <- strsplit(args$input_reduced_dim, ",")[[1]]
args$reduction_methods <- strsplit(args$reduction_methods, ",")[[1]]
args$vars_to_regress_out <- strsplit(args$vars_to_regress_out, ",")[[1]]
args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") return(TRUE)
    if (toupper(x) == "FALSE") return(FALSE)
    if (toupper(x) == "NULL") return(NULL)
  }
  return(x)
})

##  ............................................................................
##  Start                                                                   ####

sce <- read_sce(args$sce_path, read_metadata = TRUE)

sce <- reduce_dims_sce(
    sce,
    input_reduced_dim = args$input_reduced_dim,
    reduction_methods = args$reduction_methods,
    vars_to_regress_out = args$vars_to_regress_out,
    pca_dims = args$pca_dims,
    n_neighbors = args$n_neighbors,
    n_components = args$n_components,
    init = args$init,
    metric = args$metric,
    n_epochs = args$n_epochs,
    learning_rate = args$learning_rate,
    min_dist = args$min_dist,
    spread = args$spread,
    set_op_mix_ratio = args$set_op_mix_ratio,
    local_connectivity = args$local_connectivity,
    repulsion_strength = args$repulsion_strength,
    negative_sample_rate = args$negative_sample_rate,
    fast_sgd = args$fast_sgd)

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "reddim_sce"),
  write_metadata = TRUE
  )

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
