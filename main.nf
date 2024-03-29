#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/scflow
========================================================================================
    Github : https://github.com/nf-core/scflow
    Website: https://nf-co.re/scflow
    Slack  : https://nfcore.slack.com/channels/scflow
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

// Execute a single named workflow
include { SCFLOW } from './workflows/scflow'
workflow {
    SCFLOW ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
