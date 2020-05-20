#!/usr/bin/env Rscript
# Finalize SCE with manually revised celltypes
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(magrittr)

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

required$add_argument(
  "--celltype_var",
  help = "name of the column with celltype names",
  metavar = "foo/bar", 
  required = TRUE
)

required$add_argument(
  "--unique_id_var",
  help = "name of the column with unique sample ids",
  metavar = "foo/bar", 
  required = TRUE
)

required$add_argument(
  "--facet_vars",
  help = "names of variables to examine in the celltype metrics report",
  metavar = "foo/bar", 
  required = TRUE
)


required$add_argument(
  "--input_reduced_dim",
  help = "name of the reduced dimension slot to use for plots in the report",
  metavar = "foo/bar", 
  required = TRUE
)

required$add_argument(
  "--metric_vars",
  help = "names of variables to examine in the celltype metrics report",
  metavar = "foo/bar", 
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$facet_vars <- strsplit(args$facet_vars, ",")[[1]]
args$metric_vars <- strsplit(args$metric_vars, ",")[[1]]

##  ............................................................................
##  Start                                                                   ####

celltype_mappings <- read_celltype_mappings(args$celltype_mappings)

sce <- read_sce(args$sce_path)

sce <- map_custom_celltypes(
    sce, 
    mappings = celltype_mappings, 
    clusters_colname = args$clusters_colname
    )

sce <- annotate_celltype_metrics(
  sce,
  cluster_var = args$clusters_colname,
  celltype_var = args$celltype_var,
  unique_id_var = args$unique_id_var,
  facet_vars = args$facet_vars,
  input_reduced_dim = args$input_reduced_dim,
  metric_vars = args$metric_vars
)

dir.create(file.path(getwd(), "celltype_metrics_report"))

report_celltype_metrics(
  sce = sce,
  report_folder_path = file.path(getwd(), "celltype_metrics_report"),
  report_file = "scflow_celltype_metrics_report"
)

##  ............................................................................
##  Save Outputs                                                            ####

celltypes <- as.data.frame(SummarizedExperiment::colData(sce)) %>%
  dplyr::count(cluster_celltype)
colnames(celltypes) <- c("celltype", "n_cells")

write.table(
  data.frame(celltypes), 
  file = "celltypes.tsv", 
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), "final_sce")
  )

##  ............................................................................
##  Clean up                                                                ####
