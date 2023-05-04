# nf-core/scflow: Output

## Introduction

This document describes the output produced by the pipeline. Key outputs include interactive HTML reports for major analytical steps, flat-file tables, and publication-quality plots. In addition, a fully-annotated SingleCellExperiment (SCE) object is output for optional downstream tertiary analysis, in addition to individual SCEs for each sample after quality control.

The pipeline will create the directories listed below during an analysis run. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and automates a case/control scRNA-seq analysis using the following steps:

* [Check inputs](#check_inputs) - Checks the input sample sheet and manifest files (`SCFLOW_CHECKINPUTS`)
* [Quality control](#qc) - Quality control of gene-cell matrices for each individual sample (`SCFLOW_QC`)
* [Merged summary](#merged) - Quality control of merged QC tables and the merged SingleCellExperiment (SCE) object (`SCFLOW_MERGEQCTABLES` and `SCFLOW_MERGE`)
* [Integration](#integration) - Calculating latent metagene factors for the merged SCE for sample integration (`SCFLOW_INTEGRATE`, `SCFLOW_REPORTINTEGRATED`)
* [Dimension reduction](#dimension_reduction) - Dimension reduction for the merged SCE using UMAP or tSNE (`SCFLOW_REDUCEDIMS`)
* [Clustering](#clustering) - Community detection to identify clusters of cells using the Louvain/Leiden algorithm ( `SCFLOW_CLUSTER`)
* [Cell-type annotation](#celltype_annotation) - Automated cell-type annotation of the clustered SCE and identification of cell-type marker genes and calculation of relevant metrics (`SCFLOW_MAPCELLTYPES`, `SCFLOW_FINALIZE`)
* [Differential gene expression analysis](#DGE) - Performs differential gene expression analysis and generates result tables and plots (`SCFLOW_DGE`)
* [Impacted pathway analysis](#IPA) - Performs impacted pathway analysis and generates result tables and plots  (`SCFLOW_IPA`)
* [Dirichlet](#dirichlet) - Performs differential analysis of cell-type composition (`SCFLOW_DIRICHLET`)
* [Reports](#reports) - Interactive HTML reports describing results from major analytical steps of the pipeline (`SCFLOW_QC`, `SCFLOW_MERGE`, `SCFLOW_REPORTINTEGRATED`, `SCFLOW_FINALIZE`, `SCFLOW_DGE`, `SCFLOW_IPA`, `SCFLOW_DIRICHLET`)
* [Additional plots](#plots) - User-specified gene plots highlighting the expression of genes in cells plotted in reduced dimensional space
* [Additional tables](#tables) - Tables of cell-type mappings for each cluster
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution (`GET_SOFTWARE_VERSIONS`)

### Quality control

<details markdown="1">
<summary>Output files</summary>

* `quality_control/`
    * `merged.tsv` : A `.tsv` file containing detailed individual-sample QC metrics for all samples.

* `quality_cotrol/<manifest>/`
    * `qc_plot_data` : `.tsv` files for major QC values for plotting.
    * `qc_plots` : `.png` files included in the `.html` reports for major qc steps.
    * `sce/<manifest_sce>/`: Post-QC SCE for an individual sample. It is possible to use the read_sce() function of the scFlow R package to read in this object.

</details>

### Merged summary

<details markdown="1">
<summary>Output files</summary>

* `merged/`
    * `merged_plots/` : Pseudobulk plots from the `*_scflow_merged_report.html` report.
    * `merge_summary_plots/` : High-resolution plots from the `*_scflow_merged_report.html` report.

</details>

### Integration, dimension reduction and clustering

<details markdown="1">
<summary>Output files</summary>

The process `SCFLOW_REPORTINTEGRATED` saves an interactive integration and clustering HTML report to the reports/ folder.

</details>

### Celltype annotation

<details markdown="1">
<summary>Output files</summary>

* `celltype_markers/celltype_marker_plots/`
    * `*.pdf` : `.pdf` image of the marker gene plots for clusters and cluster_celltype variables.
    * `*.png` : `.png` image of the marker gene plots for clusters and cluster_celltype variables.

* `celltype_markers/celltype_marker_tables/`
    * `*.tsv` : `.tsv` files for all and top n marker genes for clusters and cluster_celltype variables

* `final/`
    * `SCE/final_sce` : The directory containing the final SCE with all metadata and dimensionality reduction, clustering and cell-type annotation. It is possible to use the read_sce() function of the scFlow R package to read in this object.
    * `celltypes.tsv` : A `.tsv` file giving the final number of all cell-types and number of nuclei/cells per cell-type.

</details>

### Differential gene expression analysis

<details markdown="1">
<summary>Output files</summary>

* `DGE/<cluster_celltype>/`
    * `*.tsv` : A `.tsv` file containing all genes with statistical results of the fitted model, including logFC, adjusted p-value, etc.
    * `de_plots` : Directory containing the volcano plot used in the `scflow_de_report.html` report.

</details>

### Impacted pathway analysis

<details markdown="1">
<summary>Output files</summary>

* `IPA/<cluster_celltype>/`
    * `<enrichment_tool>/<de_table>/` : Directory containing a `.png` plot of the top 10 significant impacted pathways and `.tsv` file containing all significantly impacted pathways.

</details>

### Dirichlet

<details markdown="1">
<summary>Output files</summary>
The process `SCFLOW_DIRICHLET` saves a differential cell-composition report to the `reports/` folder.

</details>

### Reports

<details markdown="1">
<summary>Output files</summary>

* `reports/`
    * `qc/*_scflow_qc_report.html` : Per-sample QC reports containing post-QC summaries, key parameters used, QC plots, etc.
    * `merged_report/*_scflow_merged_report.html` : The merged summary report containing inter-sample QC metrics.
    * `integration_report/integrate_report_scflow.html` : The integration report describing key parameters for integration and visual and quantitative outputs of integration performance.
    * `celltype_metrics_report/scflow_celltype_metrics_report.html` : The cell-type metrics report, including cluster and cell-type annotations, marker genes, and additional metrics.
    * `DGE/*scflow_de_report.html` : Individual reports for each differential gene expression model fit for each cell-type.
    * `IPA/*scflow_ipa_report.html` : Individual reports from impacted pathway analysis with plots and tables of  enrichment results.
    * `dirichlet_report/*dirichlet_report.html` : Dirichlet report for differential cell-type composition analysis.

</details>

### Additional plots

<details markdown="1">
<summary>Output files</summary>

* `plots/reddim_gene_plots`
    * `<plotreddim_reduction_methods>/<celltypes>` : Directories of plots of gene expression in 2D space for each gene in the `reddim_genes.yml` file.

</details>

### Additional tables

<details markdown="1">
<summary>Output files</summary>

* `tables/celltype_mappings/`
    * `celltype_mappings.tsv` : A `.tsv` file with the automated cell-type annotations generated by the process `SCFLOW_MAPCELLTYPES`.  Optionally copy this file to a new location, update it, and return to the analysis with the `--celltype_mappings` parameter to manually revise cell-type annotations for clusters.

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline and provide you with other information such as launch commands, run times, and resource usage.

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.svg`.
    * Reports generated by the pipeline: `software_versions.tsv`.

</details>
