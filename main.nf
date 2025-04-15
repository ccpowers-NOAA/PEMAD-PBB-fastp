// PROCESSES

process PRE_FASTQC {

    publishDir "${params.publishDir}/pre_fastqc", mode: 'copy'

    input:
        tuple val(sampleID), file(reads)

    output:
        path "${sampleID}_R1_fastqc.{zip,html}", emit: R1_report
        path "${sampleID}_R2_fastqc.{zip,html}", emit: R2_report

    shell:
        '''
        fastqc !{reads[0]}
        fastqc !{reads[1]}
        '''
}

process FASTP {

    publishDir "${params.publishDir}/fastp", mode: 'copy'

    input:
        tuple val(sampleID), file(reads)

    output:
<<<<<<< Updated upstream
        path "${reads[0].simpleName}_trimmed.fastq.gz", emit: R1_trimmed
        path "${reads[1].simpleName}_trimmed.fastq.gz", emit: R2_trimmed
=======
        tuple path("${reads[0].simpleName}_trimmed.fastq.gz"),
              path("${reads[1].simpleName}_trimmed.fastq.gz"), emit: trimmed_reads
>>>>>>> Stashed changes

    shell:
        '''
        fastp -i !{reads[0]} \
            -I !{reads[1]} \
            -o !{reads[0].simpleName}_trimmed.fastq.gz \
            -O !{reads[1].simpleName}_trimmed.fastq.gz
        '''
}

process POST_FASTQC {

    publishDir "${params.publishDir}/post_fastqc", mode: 'copy'

    input:
        path R1_trimmed
        path R2_trimmed

    output:
        path "${R1_trimmed.simpleName}_fastqc.{zip,html}", emit: R1_trimmed_report
        path "${R2_trimmed.simpleName}_fastqc.{zip,html}", emit: R2_trimmed_report

    shell:
        '''
        fastqc !{R1_trimmed}
        fastqc !{R2_trimmed}
        '''
}

process MULTIQC {

    publishDir "${params.publishDir}/multiqc", mode: 'copy'

    input:
        path pre_reports
        path post_reports

    output:
        path "pretrimming.html"
        path "posttrimming.html"

    shell:
        '''
        multiqc !{pre_reports} -n pretrimming.html
        multiqc !{post_reports} -n posttrimming.html
        '''
}

// WORKFLOW

workflow QC {
    // channels
    reads = Channel.fromFilePairs("${params.raw_data_dir}/SRR*_R{1,2}.fastq.gz")
    PRE_FASTQC(reads)
    FASTP(reads)
<<<<<<< Updated upstream
    POST_FASTQC(FASTP.out.R1_trimmed, FASTP.out.R2_trimmed)
=======
    POST_FASTQC(FASTP.out.trimmed_reads)
>>>>>>> Stashed changes
    MULTIQC(PRE_FASTQC.out.R1_report.mix(PRE_FASTQC.out.R2_report).collect(), 
            POST_FASTQC.out.R1_trimmed_report.mix(POST_FASTQC.out.R2_trimmed_report).collect())
}
