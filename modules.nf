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
    // echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary.csv")// from Trimming_ch


    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary1.csv"), file("${base}.scaffolds.fasta")// into Spades_Assembly_ch
    
    publishDir "${params.outdir}denovo_assembly", mode: 'copy', pattern:'*scaffolds.fasta*'

    script:

    """
    #!/bin/bash

    /SPAdes-3.15.3-Linux/bin/spades.py -s ${base}.trimmed.fastq.gz -o ${params.outdir}spades/

    cp ${base}_summary.csv ${base}_summary1.csv

    """
}
// /usr/local/bin/bbmap.sh
/*
 * STEP 2: Bam_Sorting
 * Sort bam file and collect summary statistics.
 */
process Alignment { 
    container "docker.io/buchfink/diamond:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary1.csv"), file("${base}.scaffolds.fasta")// from Spades_Assembly_ch

    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary2.csv"), file("${base}.scaffolds.fasta"), file("${base}_output.tsv"), file("${base}_all_accession.txt") //  into Alignment_ch   

    publishDir "${params.outdir}diamond_alignment", mode: 'copy', pattern:'*_output.tsv*'
    publishDir "${params.outdir}diamond_alignment", mode: 'copy', pattern:'*_all_accession.txt*'    

    script:
    """
    #!/bin/bash

    diamond blastx -d nr -q ${base}.scaffolds.fasta -o ${base}_output.tsv

    cat ${base}_output.tsv | tr "\t" "~" | cut -d"~" -f2 > ${base}_all_accession.txt

    cp ${base}_summary1.csv ${base}_summary2.csv

    """
}
/*
 * STEP 4: Analysis
 * Analysis summary creation utilizing R script.
 */
process Blast {
    container "docker.io/luciorq/entrez-direct:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
    tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary2.csv"), file("${base}.scaffolds.fasta"), file("${base}_output.tsv"), file("${base}_all_accession.txt")// from Alignment_ch    

    output:
    tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary2.csv"), file("${base}.scaffolds.fasta"), file("${base}_output.tsv"), file("${base}_all_accession.txt"), file("${base}_output.txt")//  into Blast_ch   

    publishDir "${params.outdir}blast", mode: 'copy', pattern:'*_output.txt*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_summary.csv*'

    script:
    """
    # Entrez Direct eutils Esummary of protein hits:

    cat ${base}_all_accession.txt | while read line; do esummary -db protein | xtract -pattern DocumentSummary -element Caption,TaxId, Id Title; done >> ${base}_output.txt

    cp ${base}_summary2.csv ${base}_summary.csv

    """
}
