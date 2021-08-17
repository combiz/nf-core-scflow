/*
 * Integrated pathway analysis of differentially expressed genes
 */

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SCFLOW_IPA {
    tag "${de_table_basename}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //    container 'combiz/scflow-docker:0.6.1'

    input:
    path de_table

    output:
    path '*_ipa'      , emit: ipa_results , optional: true, type: 'dir'
    path '*.html'   , emit: ipa_report  , optional: true

    script:
    de_table_basename = "${de_table.baseName}"
    def software = getSoftwareName(task.process)

    """
    export MC_CORES=${task.cpus}

    scflow_ipa.r \
    $options.args \
    --gene_file ${de_table.join(',')}

    mv ipa ${de_table_basename}_ipa

    scflow_version=\$(Rscript -e 'cat(as.character(utils::packageVersion("scFlow")))'); echo "scFlow \${scflow_version}" > "scFlow_\${scflow_version}.version.txt"
    """
}
