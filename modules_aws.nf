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
 * STEP 2: Denovo_Assembly
 * Denovo assembly of sequence reads
 */
process Denovo_Assembly {
    container "docker.io/paulrkcruz/viral_discovery_pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary.csv")// from Trimming_ch

    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary1.csv")// into Spades_Assembly_ch
    
    publishDir "${params.outdir}denovo_assembly", mode: 'copy', pattern:'*.fasta*'
    publishDir "${params.outdir}denovo_assembly", mode: 'copy', pattern:'*.info*'
    publishDir "${params.outdir}denovo_assembly", mode: 'copy', pattern:'*.log*'

    script:

    """
    #!/bin/bash

    if [ ! -d ${params.outdir}denovo_assembly ]; then
    mkdir -p ${params.outdir}denovo_assembly;
    fi;

    /root/.linuxbrew/Cellar/spades/3.15.3/bin/spades.py -t 8 -s ${base}.trimmed.fastq.gz -o ${params.outdir}denovo_assembly/ --phred-offset 33

    cp ${base}_summary.csv ${base}_summary1.csv

    """
}
// /Users/greningerlab/anaconda3/bin/spades.py
// /SPAdes-3.15.3-Linux/bin/spades.py
/*
 * STEP 3: Alignment
 * Aligns using diamond blastsx.
 */
process Alignment { 
    container "docker.io/paulrkcruz/viral_discovery_pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary1.csv")// from Spades_Assembly_ch
        file DIAMOND_DB
    output:
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary2.csv"), file("${base}_output.tsv"), file("${base}_all_accession.txt") //  into Alignment_ch   

    publishDir "${params.outdir}diamond_alignment", mode: 'copy', pattern:'*_output.tsv*'
    publishDir "${params.outdir}diamond_alignment", mode: 'copy', pattern:'*_all_accession.txt*'    

    script:
    """
    #!/bin/bash

    cp ${params.outdir}denovo_assembly/scaffolds.fasta ${params.outdir}denovo_assembly/${base}_scaffolds.fasta

    /root/.linuxbrew/Cellar/diamond/2.0.12/bin/diamond blastx -d ${DIAMOND_DB} -q ${params.outdir}denovo_assembly/${base}_scaffolds.fasta -o ${base}_output.tsv

    cat ${base}_output.tsv | tr "\t" "~" | cut -d"~" -f2 > ${base}_all_accession.txt

    cp ${base}_summary1.csv ${base}_summary2.csv

    """
}
// /opt/view/bin/diamond
/*
 * STEP 4: Generate_Summary
 * Generates the run summary using a python script.
 */
process Generate_Summary {
    container "docker.io/paulrkcruz/viral_discovery_pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input:
    tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary2.csv"), file("${base}_output.tsv"), file("${base}_all_accession.txt")// from Alignment_ch    
    file FIND_VIRUSES_PY

    output:
    tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_summary.csv"), file("${base}_output.tsv"), file("${base}_all_accession.txt"), file("${base}_output.tsv.xlsx")//  into Blast_ch   

    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*.csv*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*.xlsx*'

    script:
    """
    #!/bin/bash

    python3 ${FIND_VIRUSES_PY} ${base}_output.tsv

    cp ${base}_summary2.csv ${base}_summary.csv

    """
}