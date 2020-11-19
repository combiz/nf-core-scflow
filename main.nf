#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
========================================================================================
                         nf-core/scflow
========================================================================================
 nf-core/scflow Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/scflow
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/scflow --manifest "refs/Manifest.txt" --samplesheet "refs/SampleSheet.tsv" -params-file conf/nfx-params.json

    Mandatory arguments:
      --manifest [file]             Path to Manifest.txt file (must be surrounded with quotes)
      --samplesheet [file]          Path to SampleSheet.tsv file (must be surrounded with quotes)
      -params-file [file]           Path to nfx-params.json file
      -profile [str]                Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test, awsbatch, <institute> and more
                                    TODO: This feature will be available when configs are submitted to nf-core

    References                        If not specified in the configuration file or you wish to overwrite any of the references
      --ensembl_mappings [file]       Path to ensembl_mappings file
      --ctd_folder [path]             Path to the folder containing .ctd files for celltype annotation

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

/*
 * Create a channel for input read files
 */
 if (params.manifest) { ch_manifest = file(params.manifest, checkIfExists: true) }
 if (params.samplesheet) { ch_samplesheet = file(params.samplesheet, checkIfExists: true) }
 if (params.samplesheet) { ch_samplesheet2 = file(params.samplesheet, checkIfExists: true) } // copy for qc
 if (params.ctd_folder) { ch_ctd_folder = file(params.ctd_folder, checkIfExists: true) }
 if (params.celltype_mappings) { ch_celltype_mappings = file(params.celltype_mappings, checkIfExists: false) }
 if (params.ensembl_mappings) { ch_ensembl_mappings = file(params.ensembl_mappings, checkIfExists: false) }
 if (params.ensembl_mappings) { ch_ensembl_mappings2 = file(params.ensembl_mappings, checkIfExists: false) }
 if (params.ensembl_mappings) { ch_ensembl_mappings3 = file(params.ensembl_mappings, checkIfExists: false) }
 if (params.reddim_genes_yml) { ch_reddim_genes_yml = file(params.reddim_genes_yml, checkIfExists: false) }

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Manifest']         = params.manifest
summary['SampleSheet']      = params.samplesheet
summary['Run EmptyDrops']   = params.amb_find_cells ? "Yes" : "No"
summary['Find Singlets']    = params.mult_find_singlets ? "Yes ($params.singlets.singlets_method)" : 'No'
summary['Dimension Reds.']  = params.reddim_reduction_methods
summary['Clustering Input'] = params.clust_reduction_method
summary['DGE Method']       = params.dge_de_method == "MASTZLM" ? "$params.dge_de_method ($params.dge_mast_method)": "$params.dge_de_method"
summary['DGE Dependent Var']= params.dge_dependent_var
summary['DGE Ref Class']    = params.dge_ref_class
summary['DGE Confound Vars']= params.dge_confounding_vars
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-scflow-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/scflow Workflow Summary'
    section_href: 'https://github.com/nf-core/scflow'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    R --version > v_R.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1 - Check Inputs
 */
 process SCFLOW_CHECK_INPUTS {

  tag 'SCFLOW_CHECK_INPUTS'
  label 'process_tiny'
  //container 'google/cloud-sdk:alpine'

  echo true

  input:
    path manifest
    path samplesheet

  output:
    path 'checked_manifest.txt', emit: checked_manifest

  script:

    """
    check_inputs.r \
    --samplesheet $samplesheet \
    --manifest $manifest    
    """     

}

/*
 * STEP 2 - Single Sample QC
 */
process SCFLOW_QC {

  tag "${key}"
  label 'process_low'
  errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
  maxRetries 3
   
  echo false
   
  input:
    tuple val(key), path(mat_path)
    path samplesheet
    path ensembl_mappings

  output:
    path 'qc_report/*.html', emit: qc_report
    path 'qc_plot_data/*.tsv', emit: qc_plot_data
    path 'qc_summary/*.tsv', emit: qc_summary
    path 'qc_plots/*.png', emit: qc_plots
    path '*_sce', emit: qc_sce

  script:
    """
    scflow_qc.r \
    --samplesheet ${samplesheet} \
    --mat_path ${mat_path} \
    --key ${key} \
    --key_colname ${params.qc_key_colname} \
    --factor_vars ${params.qc_factor_vars} \
    --ensembl_mappings ${ensembl_mappings} \
    --min_library_size ${params.qc_min_library_size} \
    --max_library_size ${params.qc_max_library_size} \
    --min_features ${params.qc_min_features} \
    --max_features ${params.qc_max_features} \
    --max_mito ${params.qc_max_mito} \
    --min_ribo ${params.qc_min_ribo} \
    --max_ribo ${params.qc_max_ribo} \
    --min_counts ${params.qc_min_counts} \
    --min_cells ${params.qc_min_cells} \
    --drop_unmapped ${params.qc_drop_unmapped} \
    --drop_mito ${params.qc_drop_mito} \
    --drop_ribo ${params.qc_drop_ribo} \
    --nmads ${params.qc_nmads} \
    --find_singlets ${params.mult_find_singlets} \
    --singlets_method ${params.mult_singlets_method} \
    --vars_to_regress_out ${params.mult_vars_to_regress_out} \
    --pca_dims ${params.mult_pca_dims} \
    --var_features ${params.mult_var_features} \
    --doublet_rate ${params.mult_doublet_rate} \
    --pK ${params.mult_pK} \
    --find_cells ${params.amb_find_cells} \
    --lower ${params.amb_lower} \
    --retain ${params.amb_retain} \
    --alpha_cutoff ${params.amb_alpha_cutoff} \
    --niters ${params.amb_niters}  \
    --expect_cells ${params.amb_expect_cells}
     
    """
}

process SCFLOW_MERGE_QC_SUMMARIES {
  
  tag 'merged'
  label 'process_tiny'

  input:
    path qcs_tsv

  output:
    path '*.tsv', emit: qc_summary

  script:
    """

    merge_tables.r \
    --filepaths ${qcs_tsv.join(',')}

    """

}

process SCFLOW_MERGE {
  
  tag 'merged'
  label 'process_medium'

  input:
    path qc_passed_sces
    path ensembl_mappings

  output:
    path 'merged_sce/', emit: merged_sce
    path 'merge_plots/*.png', emit: merge_plots
    path 'merge_summary_plots/*.png', emit: merge_summary_plots    
    path 'merged_report/*.html', emit: merged_report

  script:
    """

    scflow_merge.r \
    --sce_paths ${qc_passed_sces.join(',')} \
    --ensembl_mappings ${ensembl_mappings} \
    --unique_id_var ${params.qc_key_colname} \
    --plot_vars ${params.merge_plot_vars} \
    --facet_vars ${params.merge_facet_vars} \
    --outlier_vars ${params.merge_outlier_vars}

    """

}

process SCFLOW_INTEGRATE {

  tag 'merged'
  label 'process_medium'

  input:
    path sce

  output:
    path 'integrated_sce/', emit: integrated_sce

  script:
    """

    scflow_integrate.r \
    --sce_path ${sce} \
    --method ${params.integ_method} \
    --unique_id_var ${params.integ_unique_id_var} \
    --take_gene_union ${params.integ_take_gene_union} \
    --remove_missing ${params.integ_remove_missing} \
    --num_genes ${params.integ_num_genes} \
    --combine ${params.integ_combine} \
    --keep_unique ${params.integ_keep_unique} \
    --capitalize ${params.integ_capitalize} \
    --use_cols ${params.integ_use_cols} \
    --k ${params.integ_k} \
    --lambda ${params.integ_lambda} \
    --thresh ${params.integ_thresh} \
    --max_iters ${params.integ_max_iters} \
    --nrep ${params.integ_nrep} \
    --rand_seed ${params.integ_rand_seed} \
    --knn_k ${params.integ_knn_k} \
    --k2 ${params.integ_k2} \
    --prune_thresh ${params.integ_prune_thresh} \
    --ref_dataset ${params.integ_ref_dataset} \
    --min_cells ${params.integ_min_cells} \
    --quantiles ${params.integ_quantiles} \
    --nstart ${params.integ_nstart} \
    --resolution ${params.integ_resolution} \
    --dims_use ${params.integ_dims_use} \
    --dist_use ${params.integ_dist_use} \
    --center ${params.integ_center} \
    --small_clust_thresh ${params.integ_small_clust_thresh}
    
    """

}

process SCFLOW_REDUCE_DIMS {
  
  tag 'merged'
  label 'process_medium'

  input:
    path sce

  output:
    path 'reddim_sce/', emit: reddim_sce

  script:
    """

    scflow_reduce_dims.r \
    --sce_path ${sce} \
    --input_reduced_dim ${params.reddim_input_reduced_dim} \
    --reduction_methods ${params.reddim_reduction_methods} \
    --vars_to_regress_out ${params.reddim_vars_to_regress_out} \
    --pca_dims ${params.reddim_umap_pca_dims} \
    --n_neighbors ${params.reddim_umap_n_neighbors} \
    --n_components ${params.reddim_umap_n_components} \
    --init ${params.reddim_umap_init} \
    --metric ${params.reddim_umap_metric} \
    --n_epochs ${params.reddim_umap_n_epochs} \
    --learning_rate ${params.reddim_umap_learning_rate} \
    --min_dist ${params.reddim_umap_min_dist} \
    --spread ${params.reddim_umap_spread} \
    --set_op_mix_ratio ${params.reddim_umap_set_op_mix_ratio} \
    --local_connectivity ${params.reddim_umap_local_connectivity} \
    --repulsion_strength ${params.reddim_umap_repulsion_strength} \
    --negative_sample_rate ${params.reddim_umap_negative_sample_rate} \
    --fast_sgd ${params.reddim_umap_fast_sgd} \
    --dims ${params.reddim_tsne_dims} \
    --initial_dims ${params.reddim_tsne_initial_dims} \
    --perplexity ${params.reddim_tsne_perplexity} \
    --theta ${params.reddim_tsne_theta} \
    --stop_lying_iter ${params.reddim_tsne_stop_lying_iter} \
    --mom_switch_iter ${params.reddim_tsne_mom_switch_iter} \
    --max_iter ${params.reddim_tsne_max_iter} \
    --pca_center ${params.reddim_tsne_pca_center} \
    --pca_scale ${params.reddim_tsne_pca_scale} \
    --normalize ${params.reddim_tsne_pca_normalize} \
    --momentum ${params.reddim_tsne_momentum} \
    --final_momentum ${params.reddim_tsne_final_momentum} \
    --eta ${params.reddim_tsne_eta} \
    --exaggeration_factor ${params.reddim_tsne_exaggeration_factor}
    

    """

}

process SCFLOW_CLUSTER {
  
  tag 'merged'
  label 'process_high'

  input:
    path sce

  output:
    path 'clustered_sce/', emit: clustered_sce
    path 'integration_report/', emit: integration_report

  script:
    """

    scflow_cluster.r \
    --sce_path ${sce} \
    --cluster_method ${params.clust_cluster_method} \
    --reduction_method ${params.clust_reduction_method} \
    --res ${params.clust_res} \
    --k ${params.clust_k} \
    --louvain_iter ${params.clust_louvain_iter}  \
    --categorical_covariates ${params.integ_categorical_covariates} \
    --input_reduced_dim ${params.integ_input_reduced_dim}

    """

}

process SCFLOW_MAP_CELLTYPES {
  
  tag 'merged'
  label 'process_high'

  input:
    path sce
    path ctd_folder

  output:
    path 'celltype_mapped_sce/', emit: celltype_mapped_sce
    path 'celltype_mappings.tsv', emit: celltype_mappings

  script:
    """

    scflow_map_celltypes.r \
    --sce_path ${sce} \
    --ctd_folder ${ctd_folder} \
    --clusters_colname ${params.cta_clusters_colname} \
    --cells_to_sample ${params.cta_cells_to_sample}

    """

}

process SCFLOW_FINALIZE {

  tag 'merged'
  label  'process_high'

  echo true
  
  input:
    path sce
    path celltype_mappings

  output:
    path 'final_sce/', emit: final_sce
    path 'celltypes.tsv', emit: celltypes
    path 'celltype_metrics_report', emit: celltype_metrics_report


  script:
  //def ctmappings = celltype_mappings.name != 'NO_FILE' ? "--celltype_mappings $celltype_mappings" : ''
    if( celltype_mappings.name == 'NO_FILE' )

      """
      echo "Revised celltype mappings not found: using automated celltype predictions."
      cp -R ${sce} ./final_sce
      """

    else

      """
      scflow_finalize_sce.r \
      --sce_path ${sce} \
      --celltype_mappings ${celltype_mappings} \
      --clusters_colname ${params.cta_clusters_colname} \
      --celltype_var ${params.cta_celltype_var} \
      --unique_id_var ${params.cta_unique_id_var} \
      --facet_vars ${params.cta_facet_vars.join(',')} \
      --input_reduced_dim ${params.clust_reduction_method} \
      --metric_vars ${params.cta_metric_vars}
      """

}

process SCFLOW_PLOT_REDDIM_GENES {

  label 'process_low'
   
  input:
    path sce
    path reddim_genes_yml

  output:
    path 'reddim_gene_plots/', emit: reddim_gene_plots

  script:
    
    """
    scflow_plot_reddim_genes.r \
    --sce ${sce} \
    --reduction_methods ${params.plotreddim_reduction_methods} \
    --reddim_genes_yml ${reddim_genes_yml}
     
    """
}

process SCFLOW_DGE {

  tag "${celltype} (${n_cells_str} cells) | ${de_method}"
  label 'process_medium'
  maxRetries 3
   
  input:
    path sce
    each de_method
    each ct_tuple
    path ensembl_mappings

  output:
    path 'de_table/*.tsv', emit: de_table, optional: true
    path 'de_report/*.html', emit: de_report
    path 'de_plot/*.png', emit: de_plot
    path 'de_plot_data/*.tsv', emit: de_plot_data

  script:
    celltype = ct_tuple[0]
    n_cells = ct_tuple[1].toInteger()
    n_cells_str = (Math.round(n_cells * 100) / 100000).round(1).toString() + 'k'

    """
    echo "celltype: ${celltype} n_cells: ${n_cells_str}"
    scflow_dge.r \
    --sce ${sce} \
    --celltype ${celltype} \
    --de_method ${de_method} \
    --mast_method ${params.dge_mast_method} \
    --min_counts ${params.dge_min_counts} \
    --min_cells_pc ${params.dge_min_cells_pc} \
    --rescale_numerics ${params.dge_rescale_numerics} \
    --force_run ${params.dge_force_run} \
    --pseudobulk ${params.dge_pseudobulk} \
    --celltype_var ${params.dge_celltype_var} \
    --sample_var ${params.dge_sample_var} \
    --dependent_var ${params.dge_dependent_var} \
    --ref_class ${params.dge_ref_class} \
    --confounding_vars ${params.dge_confounding_vars} \
    --random_effects_var ${params.dge_random_effects_var} \
    --fc_threshold ${params.dge_fc_threshold} \
    --ensembl_mappings ${ensembl_mappings} 
     
    """
}

process SCFLOW_IPA {

  label 'process_low'
   
  input:
    path de_table
    //each de_method
    //each celltype

  output:
    path 'ipa/**/*', optional: true, type: 'dir', emit: ipa_results
    path 'ipa/*.html', optional: true, emit: ipa_report

  script:
    """
    scflow_ipa.r \
    --gene_file ${de_table.join(',')} \
    --reference_file ${params.ipa_reference_file} \
    --enrichment_tool ${params.ipa_enrichment_tool} \
    --enrichment_method ${params.ipa_enrichment_method} \
    --enrichment_database ${params.ipa_enrichment_database}
     
    """
}


process SCFLOW_DIRICHLET {

  tag "merged"
  label  'process_low'

  echo true
  
  input:
    path sce

  output:
    path 'dirichlet_report', emit: dirichlet_report

  script:
      """
      scflow_dirichlet.r \
      --sce_path ${sce} \
      --unique_id_var ${params.dirich_unique_id_var} \
      --celltype_var ${params.dirich_celltype_var} \
      --dependent_var ${params.dirich_dependent_var} \
      --ref_class ${params.dirich_ref_class} \
      --var_order ${params.dirich_var_order} 
      """

}


workflow {  
    
  main:
    SCFLOW_CHECK_INPUTS(ch_manifest, ch_samplesheet)
    SCFLOW_QC ( SCFLOW_CHECK_INPUTS.out.checked_manifest.splitCsv(header:['key', 'filepath'], skip: 1, sep: '\t').map{ row-> tuple(row.key, row.filepath)} , ch_samplesheet2, ch_ensembl_mappings)
    //SCFLOW_QC ( CHECK_INPUTS.out.checked_manifest.splitCsv(header:true, sep: '\t').map{ row-> tuple(row.key, row.filepath)} )
    SCFLOW_MERGE_QC_SUMMARIES ( SCFLOW_QC.out.qc_summary.collect() )
    SCFLOW_MERGE ( SCFLOW_QC.out.qc_sce.collect() , ch_ensembl_mappings2 )
    SCFLOW_INTEGRATE ( SCFLOW_MERGE.out.merged_sce )
    SCFLOW_REDUCE_DIMS ( SCFLOW_INTEGRATE.out.integrated_sce )
    SCFLOW_CLUSTER ( SCFLOW_REDUCE_DIMS.out.reddim_sce )
    SCFLOW_MAP_CELLTYPES ( SCFLOW_CLUSTER.out.clustered_sce, ch_ctd_folder )
    SCFLOW_FINALIZE ( SCFLOW_MAP_CELLTYPES.out.celltype_mapped_sce, ch_celltype_mappings )
    // 
    SCFLOW_DGE( SCFLOW_FINALIZE.out.final_sce, params.dge_de_method, SCFLOW_FINALIZE.out.celltypes.splitCsv(header:['celltype', 'n_cells'], skip: 1, sep: '\t').map {row -> tuple(row.celltype, row.n_cells) } , ch_ensembl_mappings3 )
    SCFLOW_IPA( SCFLOW_DGE.out.de_table )
    //
    SCFLOW_DIRICHLET ( SCFLOW_FINALIZE.out.final_sce )
    // plotting
    SCFLOW_PLOT_REDDIM_GENES( SCFLOW_CLUSTER.out.clustered_sce, ch_reddim_genes_yml)

  /*
  publish:
    SCFLOW_CHECK_INPUTS.out.checked_manifest to: "$params.outdir/", mode: 'copy', overwrite: 'true'
    // Quality-control
    SCFLOW_QC.out.qc_report to: "$params.outdir/Reports/", mode: 'copy', overwrite: 'true'
    SCFLOW_QC.out.qc_plot_data to: "$params.outdir/Tables/Quality_Control/", mode: 'copy', overwrite: 'true'
    SCFLOW_QC.out.qc_plots to: "$params.outdir/Plots/Quality_Control/", mode: 'copy', overwrite: 'true'
    SCFLOW_QC.out.qc_sce to: "$params.outdir/SCE/Individual/", mode: 'copy', overwrite: 'true'
    SCFLOW_MERGE_QC_SUMMARIES.out.qc_summary to: "$params.outdir/Tables/Merged/", mode: 'copy', overwrite: 'true'
    // Merged SCE
    SCFLOW_MERGE.out.merged_report to: "$params.outdir/Reports/", mode: 'copy', overwrite: 'true'
    SCFLOW_MERGE.out.merge_plots to: "$params.outdir/Plots/Merged/", mode: 'copy', overwrite: 'true'
    SCFLOW_MERGE.out.merge_summary_plots to: "$params.outdir/Plots/Merged/", mode: 'copy', overwrite: 'true'
    // cluster
    SCFLOW_CLUSTER.out.integration_report to: "$params.outdir/Reports/", mode: 'copy', overwrite: 'true'
    // ct
    SCFLOW_MAP_CELLTYPES.out.celltype_mappings to: "$params.outdir/Tables/Celltype_Mappings", mode: 'copy', overwrite: 'true'
    // final
    SCFLOW_FINALIZE.out.final_sce to: "$params.outdir/SCE/", mode: 'copy', overwrite: 'true'
    SCFLOW_FINALIZE.out.celltypes to: "$params.outdir/Tables/Celltype_Mappings", mode: 'copy', overwrite: 'true'
    SCFLOW_FINALIZE.out.celltype_metrics_report to: "$params.outdir/Reports/", mode: 'copy', overwrite: 'true'
    // DE
    SCFLOW_DGE.out.de_table to: "$params.outdir/Tables/DGE", mode: 'copy', optional: true, overwrite: 'true'
    SCFLOW_DGE.out.de_report to: "$params.outdir/Reports/", mode: 'copy', overwrite: 'true'
    SCFLOW_DGE.out.de_plot to: "$params.outdir/Plots/DGE/", mode: 'copy', overwrite: 'true'
    SCFLOW_DGE.out.de_plot_data to: "$params.outdir/Tables/DGE", mode: 'copy', overwrite: 'true'
    // IPA
    SCFLOW_IPA.out.ipa_results to: "$params.outdir/Tables/", mode: 'copy', optional: true, overwrite: 'true'
    SCFLOW_IPA.out.ipa_report to: "$params.outdir/Reports/", mode: 'copy', optional: true, overwrite: 'true'
    // Dirichlet
    SCFLOW_DIRICHLET.out.dirichlet_report to: "$params.outdir/Reports/", mode: 'copy', overwrite: 'true'
    // plots
    SCFLOW_PLOT_REDDIM_GENES.out.reddim_gene_plots to: "$params.outdir/Plots/", mode: 'copy', overwrite: 'true'
*/
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/scflow] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/scflow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/scflow] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/scflow] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/scflow] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/scflow] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/scflow]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/scflow]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/scflow v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
