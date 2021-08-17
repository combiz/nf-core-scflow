#!/usr/bin/env Rscript
# Scrape versions of R package dependencies
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

# Obtain script arguments (output file path)
args <- commandArgs(trailingOnly = TRUE)
assertthat::not_empty(args)

# Retrieve package versions according with the nf-core format
pkg_versions <- tibble::tibble(
  Package = names(installed.packages()[, 3]),
  Version = paste0("v", unname(installed.packages()[, 3]))
)

# Write out package versions as a tsv file
write.table(
  pkg_versions,
  file = args[1],
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)
