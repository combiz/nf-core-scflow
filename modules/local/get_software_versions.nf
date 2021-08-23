// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    tag 'Version Info'
    label 'process_tiny'
    //cache false

    output:
    path 'software_versions.tsv'     , emit: tsv

    script: // This script is bundled with the pipeline, in nf-core/scflow/bin/
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    scrape_software_versions.r software_versions.tsv
    """
}
