/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScflow.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.manifest ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
 * Create a channel for input read files
 */
if (params.manifest) { ch_manifest = file(params.manifest, checkIfExists: true) }
if (params.input) { ch_input = file(params.input, checkIfExists: true) }
if (params.input) { ch_input2 = file(params.input, checkIfExists: true) } // copy for qc
if (params.ctd_path) { ch_ctd_path = file(params.ctd_path, checkIfExists: true) }
if (params.celltype_mappings) { ch_celltype_mappings = file(params.celltype_mappings, checkIfExists: false) }
if (params.ensembl_mappings) { ch_ensembl_mappings = file(params.ensembl_mappings, checkIfExists: false) }
if (params.ensembl_mappings) { ch_ensembl_mappings2 = file(params.ensembl_mappings, checkIfExists: false) }
if (params.ensembl_mappings) { ch_ensembl_mappings3 = file(params.ensembl_mappings, checkIfExists: false) }
if (params.reddim_genes_yml) { ch_reddim_genes_yml = file(params.reddim_genes_yml, checkIfExists: true) }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def scflow_checkinputs_options         = modules['scflow_checkinputs']
scflow_checkinputs_options.args        = ''

def scflow_qc_options                  = modules['scflow_qc']
scflow_qc_options.args                 =
    "--key_colname ${params.qc_key_colname} \
    --factor_vars ${params.qc_factor_vars} \
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
    --dpk ${params.mult_dpk} \
    --pK ${params.mult_pK} \
    --find_cells ${params.amb_find_cells} \
    --lower ${params.amb_lower} \
    --retain ${params.amb_retain} \
    --alpha_cutoff ${params.amb_alpha_cutoff} \
    --niters ${params.amb_niters}  \
    --expect_cells ${params.amb_expect_cells} \
    --species ${params.species} "

def scflow_mergeqctables_options  = modules['scflow_mergeqctables']
scflow_mergeqctables_options.args = ''

def scflow_merge_options             = modules['scflow_merge']
scflow_merge_options.args            =
    "--unique_id_var ${params.qc_key_colname} \
    --plot_vars ${params.merge_plot_vars} \
    --facet_vars ${params.merge_facet_vars} \
    --outlier_vars ${params.merge_outlier_vars} \
    --species ${params.species}"

def scflow_integrate_options         = modules['scflow_integrate']
scflow_integrate_options.args        =
    "--method ${params.integ_method} \
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
    --small_clust_thresh ${params.integ_small_clust_thresh}"

def scflow_reducedims_options        = modules['scflow_reducedims']
scflow_reducedims_options.args       =
    "--input_reduced_dim ${params.reddim_input_reduced_dim} \
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
    --normalize ${params.reddim_tsne_normalize} \
    --momentum ${params.reddim_tsne_momentum} \
    --final_momentum ${params.reddim_tsne_final_momentum} \
    --eta ${params.reddim_tsne_eta} \
    --exaggeration_factor ${params.reddim_tsne_exaggeration_factor}"

def scflow_cluster_options          = modules['scflow_cluster']
scflow_cluster_options.args         =
    "--cluster_method ${params.clust_cluster_method} \
    --reduction_method ${params.clust_reduction_method} \
    --res ${params.clust_res} \
    --k ${params.clust_k} \
    --louvain_iter ${params.clust_louvain_iter}"

def scflow_reportintegrated_options  = modules['scflow_reportintegrated']
scflow_reportintegrated_options.args =
    "--categorical_covariates ${params.integ_categorical_covariates} \
    --input_reduced_dim ${params.integ_input_reduced_dim} \
    --reddimplot_pointsize ${params.reddimplot_pointsize} \
    --reddimplot_alpha ${params.reddimplot_alpha}"

def scflow_mapcelltypes_options      = modules['scflow_mapcelltypes']
scflow_mapcelltypes_options.args     =
    "--clusters_colname ${params.cta_clusters_colname} \
    --cells_to_sample ${params.cta_cells_to_sample} \
    --species ${params.species} \
    --reddimplot_pointsize ${params.reddimplot_pointsize} \
    --reddimplot_alpha ${params.reddimplot_alpha}"

def scflow_finalize_options          = modules['scflow_finalize']
scflow_finalize_options.args         =
    "--clusters_colname ${params.cta_clusters_colname} \
    --celltype_var ${params.cta_celltype_var} \
    --unique_id_var ${params.cta_unique_id_var} \
    --facet_vars ${params.cta_facet_vars} \
    --input_reduced_dim ${params.clust_reduction_method} \
    --metric_vars ${params.cta_metric_vars} \
    --top_n ${params.cta_top_n} \
    --reddimplot_pointsize ${params.reddimplot_pointsize} \
    --reddimplot_alpha ${params.reddimplot_alpha}"

def scflow_dge_options               = modules['scflow_dge']
scflow_dge_options.args              =
    "--mast_method ${params.dge_mast_method} \
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
    --pval_cutoff ${params.dge_pval_cutoff} \
    --fc_threshold ${params.dge_fc_threshold} \
    --species ${params.species} \
    --max_cores ${params.dge_max_cores}"

def scflow_plotreddimgenes_options          = modules['scflow_plotreddimgenes']
scflow_plotreddimgenes_options.args          =
    "--reduction_methods ${params.plotreddim_reduction_methods} \
    --reddimplot_pointsize ${params.reddimplot_pointsize} \
    --reddimplot_alpha ${params.reddimplot_alpha}"

def scflow_ipa_options               = modules['scflow_ipa']
scflow_ipa_options.args              =
    "--enrichment_tool ${params.ipa_enrichment_tool} \
    --enrichment_method ${params.ipa_enrichment_method} \
    --enrichment_database ${params.ipa_enrichment_database}"

def scflow_dirichlet_options         = modules['scflow_dirichlet']
scflow_dirichlet_options.args =
    "--unique_id_var ${params.dirich_unique_id_var} \
    --celltype_var ${params.dirich_celltype_var} \
    --dependent_var ${params.dirich_dependent_var} \
    --ref_class ${params.dirich_ref_class} \
    --var_order ${params.dirich_var_order}"

def get_software_versions         = modules['get_software_versions']
get_software_versions.args = ''

include { SCFLOW_CHECKINPUTS         } from '../modules/local/process/scflow/checkinputs'       addParams( options: scflow_checkinputs_options       )
include { SCFLOW_QC                  } from '../modules/local/process/scflow/qc'                addParams( options: scflow_qc_options                )
include { SCFLOW_MERGEQCTABLES       } from '../modules/local/process/scflow/mergeqctables'     addParams( options: scflow_mergeqctables_options     )
include { SCFLOW_MERGE               } from '../modules/local/process/scflow/merge'             addParams( options: scflow_merge_options             )
include { SCFLOW_INTEGRATE           } from '../modules/local/process/scflow/integrate'         addParams( options: scflow_integrate_options         )
include { SCFLOW_REDUCEDIMS          } from '../modules/local/process/scflow/reducedims'        addParams( options: scflow_reducedims_options        )
include { SCFLOW_CLUSTER             } from '../modules/local/process/scflow/cluster'           addParams( options: scflow_cluster_options           )
include { SCFLOW_REPORTINTEGRATED    } from '../modules/local/process/scflow/reportintegrated'  addParams( options: scflow_reportintegrated_options  )
include { SCFLOW_MAPCELLTYPES        } from '../modules/local/process/scflow/mapcelltypes'      addParams( options: scflow_mapcelltypes_options      )
include { SCFLOW_FINALIZE            } from '../modules/local/process/scflow/finalize'          addParams( options: scflow_finalize_options          )
include { SCFLOW_PLOTREDDIMGENES     } from '../modules/local/process/scflow/plotreddimgenes'   addParams( options: scflow_plotreddimgenes_options   )
include { SCFLOW_DGE                 } from '../modules/local/process/scflow/dge'               addParams( options: scflow_dge_options               )
include { SCFLOW_IPA                 } from '../modules/local/process/scflow/ipa'               addParams( options: scflow_ipa_options               )
include { SCFLOW_DIRICHLET           } from '../modules/local/process/scflow/dirichlet'         addParams( options: scflow_dirichlet_options         )
include { GET_SOFTWARE_VERSIONS      } from '../modules/local/get_software_versions'            addParams( options: [publish_files : ['tsv':'']]     )


//
// MODULE: Local to the pipeline
//

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow SCFLOW {

    main:
    SCFLOW_CHECKINPUTS (
        ch_manifest,
        ch_input
    )

    SCFLOW_QC (
        SCFLOW_CHECKINPUTS.out.checked_manifest.splitCsv(
            header:['key', 'filepath'],
            skip: 1, sep: '\t'
            )
        .map { row -> tuple(row.key, row.filepath) },
        ch_input2,
        ch_ensembl_mappings
    )

    SCFLOW_MERGEQCTABLES (
        SCFLOW_QC.out.qc_summary.collect()
    )

    SCFLOW_MERGE (
        SCFLOW_QC.out.qc_sce.collect(),
        ch_ensembl_mappings2
    )

    SCFLOW_INTEGRATE (
        SCFLOW_MERGE.out.merged_sce
    )

    SCFLOW_REDUCEDIMS (
        SCFLOW_INTEGRATE.out.integrated_sce
    )

    SCFLOW_CLUSTER (
        SCFLOW_REDUCEDIMS.out.reddim_sce
    )

    SCFLOW_REPORTINTEGRATED (
        SCFLOW_CLUSTER.out.clustered_sce
    )

    SCFLOW_MAPCELLTYPES (
        SCFLOW_CLUSTER.out.clustered_sce,
        ch_ctd_path
    )

    SCFLOW_FINALIZE (
        SCFLOW_MAPCELLTYPES.out.celltype_mapped_sce,
        ch_celltype_mappings
    )

    SCFLOW_DGE (
        SCFLOW_FINALIZE.out.final_sce,
        params.dge_de_method,
        SCFLOW_FINALIZE.out.celltypes.splitCsv (
            header:['celltype', 'n_cells'], skip: 1, sep: '\t'
        )
        .map { row -> tuple(row.celltype, row.n_cells) },
        ch_ensembl_mappings3
    )

    SCFLOW_IPA (
        SCFLOW_DGE.out.de_table
    )

    SCFLOW_DIRICHLET (
        SCFLOW_FINALIZE.out.final_sce
    )

    SCFLOW_PLOTREDDIMGENES (
        SCFLOW_CLUSTER.out.clustered_sce,
        ch_reddim_genes_yml
    )

    GET_SOFTWARE_VERSIONS (
    )
}


/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
