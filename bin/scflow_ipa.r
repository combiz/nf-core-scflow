#!/usr/bin/env Rscript
# Perform impacted pathway analysis on the differential expression result table with scflow
#  Nurun Fancy <n.fancy@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = parallel::detectCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(parallel)
library(cli)

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
  "--reference_file",
  help = "full path to the reference gene file",
  metavar = ".tsv", 
  required = TRUE,
  default = "NULL"
)

required$add_argument(
  "--enrichment_tool",
  help = "name of the enrichment tool",
  metavar = "WebGestaltR,EnrichR,ROntoTools",
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
  metavar = "GO_Biological_Process,GO_Cellular_Component,GO_Molecular_Function,KEGG,Panther,Reactome,Wikipathway", 
  required = TRUE,
  default = "KEGG"
)

required$add_argument(
  "--is_output",
  help ="Whether to return output in a directory",
  metavar = "logical",
  required = TRUE,
  default = "TRUE"
)

required$add_argument(
  "--output_dir",
  help ="full path to the dir",
  metavar = "current dir",
  required = TRUE,
  default = "./"
)



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

args <- parser$parse_args()
args$enrichment_method <- strsplit(args$enrichment_method, ",")[[1]]
args$enrichment_tool <- strsplit(args$enrichment_tool, ",")[[1]]
args$enrichment_database <- strsplit(args$enrichment_database, ",")[[1]]
args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") return(TRUE)
    if (toupper(x) == "FALSE") return(FALSE)
    if (toupper(x) == "NULL") return(NULL)
  }
  return(x)
})

##  ............................................................................
##  Start impacted pathway analysis(IPA)                                                                ####

#cli::boxx(paste0("Analysing: ", args$key), float = "center")

output_dir <- file.path(getwd(), "ipa")
dir.create(output_dir)

enrichment_result <- find_impacted_pathways(gene_file = args$gene_file,
                                            enrichment_tool = args$enrichment_tool,
                                            enrichment_method = args$enrichment_method,
                                            enrichment_database = args$enrichment_database,
                                            is_output = args$is_output,
                                            output_dir = output_dir
)


report_impacted_pathway(res = enrichment_result,
                        report_folder_path = output_dir,
                        report_file = paste0(args$gene_file, "_scflow_ipa_report"))

cli::cli_text(c(
  "{cli::col_green(symbol$tick)} Analysis complete, output is found at: ", 
  "{.file {output_dir}}"
))

