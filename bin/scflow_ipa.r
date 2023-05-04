#!/usr/bin/env Rscript
# Perform impacted pathway analysis on the differential expression result
#  Nurun Fancy <n.fancy@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = parallel::detectCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(cli)
library(dplyr)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--gene_file",
  help = "full path to the gene file",
  metavar = ".tsv",
  required = TRUE,
  default = "NULL"
)

required$add_argument(
  "--enrichment_tool",
  help = "one or more enrichment tools",
  metavar = "WebGestaltR",
  required = TRUE,
  default = "WebGestaltR"
)

required$add_argument(
  "--enrichment_method",
  help = "name of the enrichment method used for webgestaltr",
  metavar = "ORA,GSEA",
  required = TRUE,
  default = "ORA"
)

required$add_argument(
  "--enrichment_database",
  help = "name of the enrichment databases",
  metavar = "GO_Biological_Process,Reactome,Wikipathway",
  required = TRUE,
  default = "GO_Biological_Process"
)

required$add_argument(
  "--species",
  help = "the biological species (e.g. mouse, human)",
  default = "human",
  required = TRUE
)

required$add_argument(
  "--logFC_threshold",
  type = "double",
  default = 0.5,
  metavar = "number",
  help = "Absolute fold-change cutoff for DE [default %(default)s]"
)

required$add_argument(
  "--padj_cutoff",
  type = "double",
  default = 0.05,
  metavar = "number",
  help = "p-value cutoff for DE [default %(default)s]"
)



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()

options("scflow_species" = args$species)

args$enrichment_method <- strsplit(args$enrichment_method, ",")[[1]]
args$enrichment_tool <- strsplit(args$enrichment_tool, ",")[[1]]
args$enrichment_database <- strsplit(args$enrichment_database, ",")[[1]]
args$gene_file <- strsplit(args$gene_file, ",")[[1]]
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

##  ............................................................................
##  Start impacted pathway analysis(IPA)                                    ####

output_dir <- file.path(getwd(), "ipa")
report_dir <- file.path(getwd())

dir.create(output_dir)
dir.create(report_dir)

for (gene_file in args$gene_file) {
  dt <- read.delim(gene_file)

  dt <- dt %>%
    dplyr::filter(
      padj <= args$padj_cutoff,
      abs(logFC) >= args$logFC_threshold
    )

  if (nrow(dt) < 5) {
    cli::cli_alert_danger("Gene list is very short!")
  } else {
    enrichment_result <- find_impacted_pathways(
      gene_file = dt,
      reference_file = NULL,
      organism = getOption("scflow_species"),
      enrichment_tool = args$enrichment_tool,
      enrichment_method = args$enrichment_method,
      enrichment_database = args$enrichment_database
    )

    for(i in names(enrichment_result)){

      output_dir_path <- paste(output_dir, i, sep = "/")
      dir.create(output_dir_path)

      res <- enrichment_result[[i]]

      lapply(
          setdiff(names(res), c("plot", "metadata")),
          function(dt) {
            write.table(res[dt],
                        file = paste(output_dir_path, "/", dt, ".tsv", sep = ""),
                        row.names = FALSE,
                        col.names = gsub(dt, "", colnames(res[[dt]])), sep = "\t")})
        lapply(
          names(res$plot),
          function(p) {
            ggplot2::ggsave(paste(output_dir_path, "/", p, ".png", sep = ""),
                            res$plot[[p]],
                            device = "png", height = 8,
                            width = 10, units = "in", dpi = 300)})
    }

    if (all(unlist(lapply(
      enrichment_result, function(dt) {
        isFALSE(dt$metadata$result)
      }
    )))) {
      cli::cli_alert_danger("No significant pathway was found at FDR 0.05")
    } else {
      report_name <- tools::file_path_sans_ext(gene_file)
      report_fp <- paste0(report_name, "_scflow_ipa_report")

      report_impacted_pathway(
        res = enrichment_result,
        report_folder_path = report_dir,
        report_file = report_fp
      )

      cli::cli_text(c(
        "{cli::col_green(symbol$tick)} Analysis complete, output is found at: ",
        "{.file {output_dir}}"
      ))
    }
  }
}
