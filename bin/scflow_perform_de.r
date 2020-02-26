#!/usr/bin/env Rscript
# Perform differential gene expression on a SingleCellExperiment Object
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = parallel::detectCores())

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
+
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
  metavar = "Microglia",
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
  type = "integer",
  default = 1,
  help = "rescale numeric variables in the model (lgl)",
  metavar = "integer",
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

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()
args$rescale_numerics <- as.logical(args$rescale_numerics)
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
  pval_cutoff = args$pval_cutoff
  )

for (result in names(de_results)) {
  write.table(de_results[[result]],
              file = paste0(
                args$celltype, "_",
                args$demethod, "_",
                result, "_DE.tsv"
              ),
              quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )
}

