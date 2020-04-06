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

// Check if genome exists in the config file
//if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
//    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
//}

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the channel below in a process, define the following:
//   input:
//   file fasta from ch_fasta
//
//params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
//if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
/*
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
*/
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */ 

 if (params.manifest) { ch_manifest = file(params.manifest, checkIfExists: true) }
 if (params.samplesheet) { ch_samplesheet = file(params.samplesheet, checkIfExists: true) }
 if (params.samplesheet) { ch_samplesheet2 = file(params.samplesheet, checkIfExists: true) } // copy for qc
 if (params.ctd_folder) { ch_ctd_folder = file(params.ctd_folder, checkIfExists: true) }
 if (params.celltype_mappings) { ch_celltype_mappings = file(params.celltype_mappings, checkIfExists: false) }
 if (params.reddim_genes_yml) { ch_reddim_genes_yml = file(params.reddim_genes_yml, checkIfExists: false) }



 /*
if (params.readPaths) {
    if (params.single_end) {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_read_files_fastqc; ch_read_files_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_read_files_fastqc; ch_read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { ch_read_files_fastqc; ch_read_files_trimming }
}
*/

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Manifest']         = params.manifest
summary['SampleSheet']      = params.samplesheet
summary['Run EmptyDrops']   = params.findcells.find_cells ? "Yes" : "No"
summary['Find Singlets']    = params.singlets.find_singlets ? "Yes ($params.singlets.singlets_method)" : 'No'
summary['Dimension Reds.']  = params.reddim.reduction_methods.join(',')
summary['Clustering Input'] = params.cluster.reduction_method
summary['DGE Method']       = params.de.de_method == "MASTZLM" ? "$params.de.de_method ($params.de.mast_method)": "$params.de.de_method"
summary['DGE Dependent Var']= params.de.dependent_var
summary['DGE Ref Class']    = params.de.ref_class.join(',')
summary['DGE Confound Vars']= params.de.confounding_vars.join(',')
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
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
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
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    //file 'software_versions_mqc.yaml' into ch_software_versions_yaml
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
 process check_inputs {

  tag "check_inputs"
  label 'process_tiny'

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
process scflow_qc {

  tag "${key}"
  label 'process_medium'
  maxRetries 3
  
  echo false
   
  input:
    tuple val(key), val(mat_path)
    path samplesheet

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
    --key_colname ${params.QC.key_colname} \
    --factor_vars ${params.QC.factor_vars.join(',')} \
    --ensembl_mappings ${params.ensembl_mappings} \
    --min_library_size ${params.QC.min_library_size} \
    --min_features ${params.QC.min_features} \
    --max_mito ${params.QC.max_mito} \
    --min_ribo ${params.QC.min_ribo} \
    --max_ribo ${params.QC.max_ribo} \
    --min_counts ${params.QC.min_counts} \
    --min_cells ${params.QC.min_cells} \
    --drop_unmapped ${params.QC.drop_unmapped} \
    --drop_mito ${params.QC.drop_mito} \
    --drop_ribo ${params.QC.drop_ribo} \
    --find_singlets ${params.singlets.find_singlets} \
    --singlets_method ${params.singlets.singlets_method} \
    --vars_to_regress_out ${params.singlets.vars_to_regress_out.join(',')} \
    --pca_dims ${params.singlets.pca_dims} \
    --var_features ${params.singlets.var_features} \
    --doublet_rate ${params.singlets.doublet_rate} \
    --pK ${params.singlets.pK} \
    --find_cells ${params.findcells.find_cells} \
    --lower ${params.findcells.lower} \
    --retain ${params.findcells.retain} \
    --alpha_cutoff ${params.findcells.alpha_cutoff} \
    --niters ${params.findcells.niters}
     
    """
}

process merge_qc_summaries {
  
  tag "merged"
  label 'process_tiny'

  input:
    path( qcs_tsv )

  output:
    path '*.tsv', emit: qc_summary

  script:
    """

    merge_tables.r \
    --filepaths ${qcs_tsv.join(',')}

    """

}

process scflow_merge {
  
  tag "merged"
  label 'process_medium'

  input:
    path( qc_passed_sces )

  output:
    path 'merged_sce/', emit: merged_sce
    path 'merge_plots/*.png', emit: merge_plots
    path 'merge_summary_plots/*.png', emit: merge_summary_plots    
    path 'merged_report/*.html', emit: merged_report

  script:
    """

    scflow_merge.r \
    --sce_paths ${qc_passed_sces.join(',')} \
    --ensembl_mappings ${params.ensembl_mappings} \
    --unique_id_var ${params.QC.key_colname} \
    --plot_vars ${params.merge.plot_vars.join(',')} \
    --facet_vars ${params.merge.facet_vars.join(',')} \
    --outlier_vars ${params.merge.outlier_vars.join(',')} \
    --outlier_mads ${params.merge.outlier_mads}

    """

}

process scflow_integrate {

  tag "merged"
  label 'process_medium'

  input:
    path( sce )

  output:
    path 'integrated_sce/', emit: integrated_sce

  script:
    """

    scflow_integrate.r \
    --sce_path ${sce} \
    --method ${params.integrate.method} \
    --unique_id_var ${params.integrate.unique_id_var} \
    --take_gene_union ${params.integrate.take_gene_union} \
    --remove_missing ${params.integrate.remove_missing} \
    --num_genes ${params.integrate.num_genes} \
    --combine ${params.integrate.combine} \
    --keep_unique ${params.integrate.keep_unique} \
    --capitalize ${params.integrate.capitalize} \
    --do_plot ${params.integrate.do_plot} \
    --cex_use ${params.integrate.cex_use} \
    --use_cols ${params.integrate.use_cols} \
    --k ${params.integrate.k} \
    --lambda ${params.integrate.lambda} \
    --thresh ${params.integrate.thresh} \
    --max_iters ${params.integrate.max_iters} \
    --nrep ${params.integrate.nrep} \
    --h_init ${params.integrate.h_init} \
    --w_init ${params.integrate.w_init} \
	  --v_init ${params.integrate.v_init} \
    --rand_seed ${params.integrate.rand_seed} \
    --print_obj ${params.integrate.print_obj} \
    --knn_k ${params.integrate.knn_k} \
    --k2 ${params.integrate.k2} \
    --prune_thresh ${params.integrate.prune_thresh} \
    --ref_dataset ${params.integrate.ref_dataset} \
    --min_cells ${params.integrate.min_cells} \
    --quantiles ${params.integrate.quantiles} \
    --nstart ${params.integrate.nstart} \
    --resolution ${params.integrate.resolution} \
	  --dims_use ${params.integrate.dims_use} \
    --dist_use ${params.integrate.dist_use} \
    --center ${params.integrate.center} \
    --small_clust_thresh ${params.integrate.small_clust_thresh} \
    --id_number ${params.integrate.id_number} \
    --print_mod ${params.integrate.print_mod} \
    --print_align_summary ${params.integrate.print_align_summary}
	
	"""

}

process scflow_reduce_dims {
  
  tag "merged"
  label 'process_low'

  input:
    path( sce )

  output:
    path 'reddim_sce/', emit: reddim_sce

  script:
    """

    scflow_reduce_dims.r \
    --sce_path ${sce} \
    --input_reduced_dim ${params.reddim.input_reduced_dim.join(',')} \
    --reduction_methods ${params.reddim.reduction_methods.join(',')} \
    --vars_to_regress_out ${params.reddim.vars_to_regress_out.join(',')} \
    --pca_dims ${params.reddim.pca_dims} \
    --n_neighbors ${params.reddim.n_neighbors} \
    --n_components ${params.reddim.n_components} \
    --init ${params.reddim.init} \
    --metric ${params.reddim.metric} \
    --n_epochs ${params.reddim.n_epochs} \
    --learning_rate ${params.reddim.learning_rate} \
    --min_dist ${params.reddim.min_dist} \
    --spread ${params.reddim.spread} \
    --set_op_mix_ratio ${params.reddim.set_op_mix_ratio} \
    --local_connectivity ${params.reddim.local_connectivity} \
    --repulsion_strength ${params.reddim.repulsion_strength} \
    --negative_sample_rate ${params.reddim.negative_sample_rate} \
    --fast_sgd ${params.reddim.fast_sgd}

    """

}

process scflow_cluster {
  
  tag "merged"
  label 'process_medium'

  input:
    path( sce )

  output:
    path 'clustered_sce/', emit: clustered_sce

  script:
    """

    scflow_cluster.r \
    --sce_path ${sce} \
    --cluster_method ${params.cluster.cluster_method} \
    --reduction_method ${params.cluster.reduction_method} \
    --res ${params.cluster.res} \
    --k ${params.cluster.k} \
    --louvain_iter ${params.cluster.louvain_iter}

    """

}

process scflow_map_celltypes {
  
  tag "merged"
  label 'process_medium'

  input:
    path( sce )
    path ctd_folder

  output:
    path 'celltype_mapped_sce/', emit: celltype_mapped_sce
    path 'celltype_mappings.tsv', emit: celltype_mappings

  script:
    """

    scflow_map_celltypes.r \
    --sce_path ${sce} \
    --ctd_folder ${ctd_folder} \
    --clusters_colname ${params.mapct.clusters_colname} \
    --cells_to_sample ${params.mapct.cells_to_sample}

    """

}

process scflow_finalize {

  tag "merged"
  label  'process_low'

  echo true
  
  input:
    path (sce)
    path celltype_mappings

  output:
    path 'final_sce/', emit: final_sce
    path 'celltypes.tsv', emit: celltypes


  script:
  //def ctmappings = celltype_mappings.name != 'NO_FILE' ? "--celltype_mappings $celltype_mappings" : ''
    if( celltype_mappings.name == 'NO_FILE' )

      """
      echo "Revised celltype mappings not found: using automated celltype predictions."
      cp -R ${sce} ./final_sce
      """

    else

      """
      echo "Applying revised celltype mappings.."
      scflow_finalize_sce.r \
      --sce_path ${sce} \
      --celltype_mappings ${celltype_mappings} \
      --clusters_colname ${params.mapct.clusters_colname}
      """

}

process scflow_plot_reddim_genes {

  label 'process_low'
   
  input:
    path( sce )
    path ( reddim_genes_yml )

  output:
    path 'reddim_gene_plots/', emit: reddim_gene_plots

  script:
    
    """
    scflow_plot_reddim_genes.r \
    --sce ${sce} \
    --reduction_methods ${params.plot_reddim_genes.reduction_methods.join(',')} \
    --reddim_genes_yml ${reddim_genes_yml}
     
    """
}

process scflow_perform_de {

  tag "${celltype} (${n_cells_str} cells) | ${de_method}"
  label 'process_medium'
  maxRetries 3
  
  echo true
   
  input:
    path( sce )
    each de_method
    each ct_tuple

  output:
    path '*.tsv', emit: de_table, optional: true

  script:
    celltype = ct_tuple[0]
    n_cells = ct_tuple[1].toInteger()
    n_cells_str = (Math.round(n_cells * 100) / 100000).round(1).toString() + 'k'

    """
    echo "celltype: ${celltype} n_cells: ${n_cells_str}"
    scflow_perform_de.r \
    --sce ${sce} \
    --celltype ${celltype} \
    --de_method ${de_method} \
    --mast_method ${params.de.mast_method} \
    --min_counts ${params.de.min_counts} \
    --min_cells_pc ${params.de.min_cells_pc} \
    --rescale_numerics ${params.de.rescale_numerics} \
    --force_run ${params.de.force_run} \
    --dependent_var ${params.de.dependent_var} \
    --ref_class ${params.de.ref_class} \
    --confounding_vars ${params.de.confounding_vars.join(',')} \
    --random_effects_var ${params.de.random_effects_var} \
    --fc_threshold ${params.de.fc_threshold}
     
    """
}

process scflow_perform_ipa {

  label 'process_medium'
   
  input:
    path( de_table )
    //each de_method
    //each celltype

  output:
    path 'ipa/**/*', optional: true, type: 'dir', emit: ipa_results
    path 'ipa/*.html', optional: true, emit: ipa_report

  script:
    """
    scflow_ipa.r \
    --gene_file ${de_table.join(',')} \
    --reference_file ${params.IPA.reference_file} \
    --enrichment_tool ${params.IPA.enrichment_tool.join(',')} \
    --enrichment_method ${params.IPA.enrichment_method} \
    --enrichment_database ${params.IPA.enrichment_database.join(',')} \
    --is_output ${params.IPA.is_output} \
    --output_dir ${params.IPA.output_dir}
     
    """
}

process scflow_traject {

  label 'process_low'
  
  echo true
   
  input:
    path( sce )

  output:
    //path '*.tsv', emit: de_table

  script:
    """
    echo hello world
     
    """
}

workflow {  
    
  main:
    check_inputs(ch_manifest, ch_samplesheet)
    scflow_qc ( check_inputs.out.checked_manifest.splitCsv(header:['key', 'filepath'], skip: 1, sep: '\t').map{ row-> tuple(row.key, row.filepath)} , ch_samplesheet2 )
    //scflow_qc ( check_inputs.out.checked_manifest.splitCsv(header:true, sep: '\t').map{ row-> tuple(row.key, row.filepath)} )
    merge_qc_summaries ( scflow_qc.out.qc_summary.collect() )
    scflow_merge ( scflow_qc.out.qc_sce.collect() )
	  scflow_integrate ( scflow_merge.out.merged_sce )
    scflow_reduce_dims ( scflow_integrate.out.integrated_sce )
    scflow_cluster ( scflow_reduce_dims.out.reddim_sce )
    scflow_map_celltypes ( scflow_cluster.out.clustered_sce, ch_ctd_folder )
    scflow_finalize ( scflow_map_celltypes.out.celltype_mapped_sce, ch_celltype_mappings )
    // 
    scflow_perform_de( scflow_finalize.out.final_sce, params.de.de_method, scflow_finalize.out.celltypes.splitCsv(header:['celltype', 'n_cells'], skip: 1, sep: '\t').map {row -> tuple(row.celltype, row.n_cells) } )
    scflow_perform_ipa( scflow_perform_de.out.de_table )
    //
    scflow_traject( scflow_finalize.out.final_sce )
    // plotting
    scflow_plot_reddim_genes( scflow_finalize.out.final_sce, ch_reddim_genes_yml)

  
  publish:
    check_inputs.out.checked_manifest to: "$params.outdir/", mode: 'copy'
    // Quality-control
    scflow_qc.out.qc_report to: "$params.outdir/qc/", mode: 'copy'
    scflow_qc.out.qc_plot_data to: "$params.outdir/qc/", mode: 'copy'
    scflow_qc.out.qc_plots to: "$params.outdir/qc/", mode: 'copy'
    scflow_qc.out.qc_sce to: "$params.outdir/sce", mode: 'copy'
    merge_qc_summaries.out.qc_summary to: "$params.outdir/qc", mode: 'copy'
    // Merged SCE
    scflow_merge.out.merged_report to: "$params.outdir/merge/", mode: 'copy'
    scflow_merge.out.merge_plots to: "$params.outdir/merge/plots", mode: 'copy'
    scflow_merge.out.merge_summary_plots to: "$params.outdir/merge/plots/summary_plots", mode: 'copy'
    // ct
    scflow_map_celltypes.out.celltype_mappings to: "$params.outdir/celltype_mappings", mode: 'copy'
    // final
    scflow_finalize.out.final_sce to: "$params.outdir/", mode: 'copy', overwrite: 'true'
    scflow_finalize.out.celltypes to: "$params.outdir/celltype_mappings", mode: 'copy', overwrite: 'true'
    // DE
    scflow_perform_de.out.de_table to: "$params.outdir/de", mode: 'copy', optional: true
    // IPA
    scflow_perform_ipa.out.ipa_results to: "$params.outdir/", mode: 'copy', optional: true
    scflow_perform_ipa.out.ipa_report to: "$params.outdir/", mode: 'copy', optional: true
    // plots
    scflow_plot_reddim_genes.out.reddim_gene_plots to: "$params.outdir/", mode: 'copy'

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
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir"]
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
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
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
