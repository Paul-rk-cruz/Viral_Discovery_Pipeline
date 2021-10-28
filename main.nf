#!/usr/bin/env nextflow

/*
========================================================================================
                        Viral Discovery Pipeline v1.0
========================================================================================
Github Repo:
Greninger Lab
https://github.com/Paul-rk-cruz/Viral_Discovery_Pipeline

Author:
Paul RK Cruz <kurtisc@uw.edu>
Kate Juergens <katej16@uw.edu>

UW Medicine | Virology
Department of Laboratory Medicine and Pathology
University of Washington
Created: October 12, 2021
Updated: October 12, 2021
LICENSE: GNU
----------------------------------------------------------------------------------------
*/
// Nextflow dsl2
nextflow.enable.dsl=2
// Set pipeline version
version = '1.0'
params.helpMsg = false
def helpMsg() {
    log.info"""
	 _______________________________________________________________________________
     Viral Discovery Pipeline :  Version ${version}
	________________________________________________________________________________
    
	Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run FILE_PATH/Viral_Discovery_Pipeline/main.nf --input PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR --SingleEnd
    Valid CLI Arguments:
    REQUIRED:
      --input                       Path to input fastq.gz folder
      --outdir                      The output directory where the results will be saved

    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
// Setup Parameters to default values
params.input = false
params.outdir = false
params.SETTING = "2:30:10:1:true"
params.LEADING = "3"
params.TRAILING = "3"
params.SWINDOW = "4:30"
params.MINLEN = "75"
// Check if input is set
if (params.input == false) {
    println( "Must provide an input directory with --input") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.input.endsWith("/")){
   params.input = "${params.input}/"
}
// if OUTDIR not set
if (params.outdir == false) {
    println( "Must provide an output directory with --outdir") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.outdir.endsWith("/")){
   params.outdir = "${params.outdir}/"
}
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException("> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version $nextflow_req_v")
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}
// Setup file paths
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
// Workflow display header
def hpvheader() {
    return """
    """.stripIndent()
}
// log files header
// log.info hpvheader()
log.info "_______________________________________________________________________________"
log.info " Viral Discovery Pipeline :  v${version}"
log.info "_______________________________________________________________________________"
def summary = [:]
summary['Configuration Profile:'] = workflow.profile
summary['Current directory path:']        = "$PWD"
summary['Pipeline directory path:']          = workflow.projectDir
summary['Input directory path:']               = params.input
summary['Output directory path:']          = params.outdir
summary['Work directory path:']         = workflow.workDir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Trimmomatic adapters:'] = ADAPTERS_SE
summary["Trimmomatic read length (minimum):"] = params.MINLEN
summary["Trimmomatic Setting:"] = params.SETTING
summary["Trimmomatic Sliding Window:"] = params.SWINDOW
summary["Trimmomatic Leading:"] = params.LEADING
summary["Trimmomatic Trailing:"] = params.TRAILING
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "_______________________________________________________________________________"


//
// Import processes
// 

include { Trimming } from './modules.nf'
include { Denovo_Assembly } from './modules.nf'
include { Alignment } from './modules.nf'
include { Blast } from './modules.nf'

// Create channel for input reads: single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.input}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.input} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.input}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    Trimming (
        input_read_ch,
        ADAPTERS_SE,
        params.MINLEN,
        params.SETTING,
        params.LEADING,
        params.TRAILING,
        params.SWINDOW
    )
    Denovo_Assembly (
        Trimming.out[0],
    )
    Alignment (
        Denovo_Assembly.out[0]
    )
    Blast (
        Alignment.out[0]
    )
}
