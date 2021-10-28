// HPV Pipeline Workflow Processes
/*
 * STEP 1: Trim_Reads
 * Trimming of low quality and short NGS sequences.
 */
process Trimming {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        file R1 //from input_read_ch
        file ADAPTERS_SE
        val MINLEN
        val SETTING
        val LEADING
        val TRAILING
        val SWINDOW
    output:
        tuple env(base),file("*.trimmed.fastq.gz"), file("*summary.csv")// into Trimming_ch

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    
    base=`basename ${R1} ".fastq.gz"`
    echo \$base
    
    /usr/local/miniconda/bin/trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS_SE}:${SETTING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SWINDOW} MINLEN:${MINLEN}
    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    printf "\$num_trimmed" >> ${R1}_num_trimmed.txt
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed> \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    ls -latr

    """
}
/*
 * STEP 2: Alignment
 * Align NGS Sequence reads to HPV-ALL multifasta
 */
process Denovo_Assembly {
    container "docker.io/staphb/spades:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary.csv")// from Trimming_ch


    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_stats1.csv")// into Align_ch
    
    publishDir "${params.outdir}bbmap_scaf_stats", mode: 'copy', pattern:'*_hpvAll_scafstats.txt*'
 

    script:

    """
    #!/bin/bash

    spades.py -s input file -o output directory

    cp ${base}_summary.csv ${base}_stats1.csv

    """
}
// /usr/local/bin/bbmap.sh
/*
 * STEP 2: Bam_Sorting
 * Sort bam file and collect summary statistics.
 */
process Alignment { 
    container "quay.io/greninger-lab/swift-pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz")// from Align_ch

    output:
        tuple val(base), file("${base}_trim_stats.csv")//  into Analysis_ch   

    publishDir "${params.outdir}bam_sorted", mode: 'copy', pattern:'*_hpvAll.sorted.bam*'

    script:
    """
    #!/bin/bash

    diamond blastx -d nr -q contigs.fasta -o output.tsv



    cp ${base}_stats1.csv ${base}_trim_stats.csv

    """
}
/*
 * STEP 4: Analysis
 * Analysis summary creation utilizing R script.
 */
process Analysis {
    container "docker.io/rocker/tidyverse:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
    file("${base}_hpvAll_scafstats.txt")// from Bbmap_scaf_stats_ch.collect()     

    script:
    """


    """
}