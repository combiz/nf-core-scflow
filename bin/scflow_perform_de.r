#!/usr/bin/env Rscript
# Perform differential gene expression on a SingleCellExperiment Object
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

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
  "--sce",
  help = "path to SingleCellExperiment directory",
  metavar = "/dir/sce/",
  required = TRUE
)

required$add_argument(
  "--celltype",
  help = "celltype to subset for DE analysis",
  metavar = "DESEQ2PB",
  required = TRUE
)

required$add_argument(
  "--de_method",
  help = "differential gene expression method",
  metavar = "MAST",
  required = TRUE
)

required$add_argument(
  "--mast_method",
  help = "differential gene expression sub-method for MAST",
  metavar = "bayesglm",
  required = TRUE
)

required$add_argument(
  "--min_counts",
  type = "integer",
  default = 1,
  help = "minimum library size (counts) per cell",
  metavar = "N",
  required = TRUE
)

required$add_argument(
  "--min_cells_pc",
  type = "double",
  default = 0.10,
  metavar = "N",
  help = "minimum percentage of cells with min_counts"
)

required$add_argument(
  "--rescale_numerics",
  help = "rescale numeric variables in the model (lgl)",
  metavar = "TRUE",
  required = TRUE
)

required$add_argument(
  "--pseudobulk",
  help = "perform pseudobulking option (lgl)",
  metavar = "TRUE",
  required = TRUE
)

required$add_argument(
  "--celltype_var",
  help = "celltype variable",
  metavar = "cluster_celltype",
  required = TRUE
)

required$add_argument(
  "--sample_var",
  help = "sample variable",
  metavar = "manifest",
  required = TRUE
)

required$add_argument(
  "--force_run",
  help = "force run if non-full-rank (lgl)",
  metavar = "TRUE",
  required = TRUE
)

required$add_argument(
  "--dependent_var",
  help = "dependent variable",
  metavar = "group",
  required = TRUE
)

required$add_argument(
  "--ref_class",
  help = "reference class within dependent variable",
  metavar = "Control",
  required = TRUE
)

required$add_argument(
  "--confounding_vars",
  help = "confounding variables",
  metavar = "age,sex,pc_mito",
  required = TRUE
)

required$add_argument(
  "--random_effects_var",
  help = "random effects variable",
  metavar = "individual",
  required = TRUE
)

required$add_argument(
  "--fc_threshold",
  type = "double",
  default = 1.1,
  metavar = "number",
  help = "Absolute fold-change cutoff for DE [default %(default)s]"
)

required$add_argument(
  "--pval_cutoff",
  type = "double",
  default = 0.05,
  metavar = "number",
  help = "p-value cutoff for DE [default %(default)s]"
)

required$add_argument(
  "--ensembl_mappings",
  help = "path to ensembl mappings file",
  metavar = "tsv", 
  required = TRUE
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()
args$rescale_numerics <- as.logical(args$rescale_numerics)
args$pseudobulk <- as.logical(args$pseudobulk)
args$force_run <- as.logical(args$force_run)
if(tolower(args$random_effects_var) == "null") args$random_effects_var <- NULL
args$confounding_vars <- strsplit(args$confounding_vars, ",")[[1]]

#   ____________________________________________________________________________
#   Start DE                                                                ####

write(sprintf(
  "##### Starting DE of %s cells with %s",
  args$celltype, args$demethod
), stdout())

sce <- read_sce(args$sce)

sce_subset <- sce[, sce$cluster_celltype == args$celltype]

if (args$pseudobulk) {
  pb_str <- pb_str <- "_pb"
  sce_subset <- pseudobulk_sce(
    sce_subset,
    keep_vars = c(args$dependent_var, args$confounding_vars, args$random_effects_var),
    assay_name = "counts",
    celltype_var = args$celltype_var,
    sample_var = args$sample_var
  )
} else {
  pb_str <- ""
}

de_results <- perform_de(
  sce_subset,
  de_method = args$de_method,
  min_cells_pc = args$min_cells_pc,
  rescale_numerics = args$rescale_numerics,
  dependent_var = args$dependent_var,
  ref_class = args$ref_class,
  confounding_vars = args$confounding_vars,
  random_effects_var = args$random_effects_var,
  fc_threshold = args$fc_threshold,
  pval_cutoff = args$pval_cutoff,
  mast_method = args$mast_method,
  force_run = args$force_run,
  ensembl_mapping_file = args$ensembl_mappings
  )

for (result in names(de_results)) {
  if (dim(de_results[[result]])[[1]] > 0) {
    write.table(de_results[[result]],
                file = paste0(
                  args$celltype, "_",
                  args$demethod, "_",
                  pb_str,
                  result, "_DE.tsv"
                ),
                quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
    )
  } else { 
    print(sprintf("No DE genes found for %s", result)) 
  }  
}


