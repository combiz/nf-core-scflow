#!/usr/bin/env Rscript
# Perform quality-control on a feature-barcode matrix with scflow
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

#options(mc.cores = parallel::detectCores())

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
  "--samplesheet",
  help = "full path to the sample sheet tsv file",
  metavar = "SampleSheet.tsv", 
  required = TRUE
)

required$add_argument(
  "--key_colname",
  help = "sample sheet column name with unique sample identifiers",
  metavar = "manifest", 
  required = TRUE
)

required$add_argument(
  "--key",
  help = "unique identifier in sample sheet column specified by key_colname",
  metavar = "hirol", 
  required = TRUE
)

required$add_argument(
  "--mat_path",
  help = "folder path of sparse matrix (cellranger output)",
  metavar = "out", 
  required = TRUE
)

required$add_argument(
  "--ensembl_mappings",
  help = "path to ensembl mappings file",
  metavar = "tsv", 
  required = TRUE
)

required$add_argument(
  "--min_library_size",
  type = "integer", 
  default = 300,
  help = "minimum library size (counts) per cell",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--min_features",
  type = "integer", 
  default = 100,
  help = "minimum features (expressive genes) per cell",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--max_mito",
  type = "double", 
  default = 0.1,
  help = "maximum proportion of counts mapping to mitochondrial genes",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--min_ribo",
  type = "double", 
  default = 0.0,
  help = "minimum proportion of counts mapping to ribosomal genes",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--max_ribo",
  type = "double", 
  default = 1.0,
  help = "maximum proportion of counts mapping to ribosomal genes",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--min_counts",
  type = "integer", 
  default = 2,
  help = "expressive genes must have >=min_counts in >=min_cells",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--min_cells",
  type = "integer", 
  default = 2,
  help = "expressive genes must have >=min_counts in >=min_cells",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--drop_unmapped",
  type = "integer", 
  default = 1,
  help = "drop genes which could not be mapped to gene names (lgl)",
  metavar = "integer", 
  required = TRUE
)

required$add_argument(
  "--drop_mito",
  type = "integer", 
  default = 1,
  help = "drop mitochondrial genes (lgl)",
  metavar = "integer", 
  required = TRUE
)

required$add_argument(
  "--drop_ribo",
  type = "integer", 
  default = 0,
  help = "drop ribosomal genes (lgl)",
  metavar = "integer", 
  required = TRUE
)

required$add_argument(
  "--find_singlets",
  type = "integer", 
  default = 1,
  help = "run a singlet finding algorithm (lgl)",
  metavar = "integer", 
  required = TRUE
)

required$add_argument(
  "--singlets_method",
  help = "method to identify singlets",
  metavar = "doubletfinder", 
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
  default = 10,
  help = "number of principal components for singlet finding",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--var_features",
  type = "integer", 
  default = 2000,
  help = "number of variable features for singlet finding",
  metavar = "N", 
  required = TRUE
)

required$add_argument(
  "--doublet_rate",
  type = "double", 
  default = 0.075,
  help = "estimated doublet rate",
  metavar = "N", 
  required = TRUE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args[startsWith(names(args), "drop_")] <- 
  as.logical(args[startsWith(names(args), "drop_")])
args$find_singlets <- as.logical(args$find_singlets)
args$vars_to_regress_out <- strsplit(args$vars_to_regress_out, ",")[[1]]

##  ............................................................................
##  Start QC                                                                ####

cli::boxx(paste0("Analysing: ", args$key), float = "center")

mat <- scflow::read_sparse_matrix(args$mat_path)

metadata <- read_metadata(
  unique_key = args$key,
  key_colname = args$key_colname,
  samplesheet_path = args$samplesheet
)

sce <- generate_sce(mat, metadata)

rm(mat, metadata)

sce <- annotate_sce(
  sce = sce,
  min_library_size = args$min_library_size,
  min_features = args$min_features,
  max_mito = args$max_mito,
  min_ribo = args$min_ribo,
  max_ribo = args$max_ribo,
  min_counts = args$min_counts,
  min_cells = args$min_cells,
  drop_unmapped = args$drop_unmapped,
  drop_mito = args$drop_mito,
  drop_ribo = args$drop_ribo,
  annotate_genes = TRUE,
  annotate_cells = TRUE,
  ensembl_mapping_file = args$ensembl_mappings
)

sce <- filter_sce(
  sce, 
  filter_genes = TRUE, 
  filter_cells = TRUE
)

if(args$find_singlets) {
  sce <- find_singlets(
    sce = sce,
    singlet_find_method = args$singlets_method,
    vars_to_regress_out = args$vars_to_regress_out,
    pca_dims = args$pca_dims,
    var_features = args$var_features,
    doublet_rate = args$doublet_rate
  )

  sce <- filter_sce(
    sce, 
    filter_genes = TRUE, 
    filter_cells = TRUE
  )
}

dir.create(file.path(getwd(), "qc_report"))

report_qc_sce(
  sce = sce,
  report_folder_path = file.path(getwd(), "qc_report"),
  report_file = paste0(args$key, "_scflow_qc_report")
)

print("Analysis complete, saving outputs..")

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), paste0(args$key, "_sce"))
  )

new_dirs <- c(
  "qc_plot_data",
  "qc_summary",
  "qc_plots")

#make dirs
purrr::walk(new_dirs, ~ dir.create(file.path(getwd(), .)))

# Save QC plots (tables)
for (df in names(sce@metadata$qc_plot_data)) {
  write.table(
    sce@metadata$qc_plot_data[[df]],
    file.path(getwd(), "qc_plot_data", 
              paste0(args$key, "_", df, ".tsv")),
    sep = "\t",
    col.names = TRUE, row.names = FALSE)
}

# Save QC summary table
write.table(
  cbind(sce@metadata$metadata, sce@metadata$qc_summary),
  file.path(getwd(), "qc_summary", 
            paste0(args$key, "_qc_summary.tsv")),
  sep = "\t",
  col.names = TRUE, row.names = FALSE)

# Save QC plots (images)
for(pname in names(sce@metadata$qc_plots)) {
  png(file.path(getwd(), "qc_plots", 
                paste0(args$key, "_", pname, ".png")), 
      width = 247, height = 170, units = "mm", res = 600)
  print(sce@metadata$qc_plots[[pname]])
  dev.off()
}

# Save doublet finder plots, square
for(pname in names(sce@metadata$qc_plots$doublet_finder)) {
  png(file.path(getwd(), "qc_plots", 
                paste0(args$key, "_", pname, "_doublet_finder.png")), 
      width = 170, height = 170, units = "mm", res = 600)
  print(sce@metadata$qc_plots$doublet_finder[[pname]])
  dev.off()
}

##  ............................................................................
##  Clean up                                                                ####

# Clear biomart cache
biomaRt::biomartCacheClear()
