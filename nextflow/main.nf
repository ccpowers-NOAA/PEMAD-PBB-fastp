// PROCESSES

process PRE_FASTQC {
    executor params.pre_fastqc.executor

    publishDir "${launchDir}/${params.publishDir}/pre_fastqc", mode: 'copy'

    input:
        tuple val(sampleID), file(reads)

    output:
        path "${sampleID}_R1_fastqc.zip"
        path "${sampleID}_R1_fastqc.html"
        path "${sampleID}_R2_fastqc.zip"
        path "${sampleID}_R2_fastqz.html"

    shell:
        '''
        fastqc ${reads[0]}
        fastqc ${reads[1]}
        '''
}

process FASTP {

    publishDir "${launchDir}/${params.publishDir}/fastp", mode: 'copy'

    input:
        path input_reads

    output:
        path "${reads[0].simpleName}_trimmed.fastq.gz", emit: R1_trimmed
        path "${reads[1].simpleName}_trimmed.fastq.gz", emit: R2_trimmed

    shell:
        '''
        fastp -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${reads[0].simpleName}_trimmed.fastq.gz \
            -O ${reads[1].simpleName}_trimmed.fastq.gz
        '''
}

process POST_FASTQC {

    publishDir "${launchDir}/${params.publishDir}/post_fastqc", mode: 'copy'

    input:
        path R1_trimmed
        path R2_trimmed

    output:
        path "${R1_trimmed.simpleName}_fastqc.zip"
        path "${R1_trimmed.simpleName}_fastqc.html"
        path "${R2_trimmed.simpleName}_fastqc.zip"
        path "${R2_trimmed.simpleName}_fastqc.html"

    shell:
        '''
        fastqc ${R1_trimmed}
        fastqc ${R2_trimmed}
        '''
}

process MULTIQC {

    publishDir "${launchDir}/${params.publishDir}/multiqc", mode: 'copy'

    output:
        path "pretrimming.html"
        path "posttrimming.html"

    shell:
        '''
        multiqc ${launchDir}/${params.publishDir}/pre_fastqc -n pretrimming.html
        multiqc ${launchDir}/${params.publishDir}/post_fastqc -n posttrimming.html
        '''
}


// WORKFLOW

workflow QC {
    // channels
    reads = Channel.fromFilePairs("${params.raw_data_dir}/SRR*_R{1,2}.fastq.gz")
    PRE_FASTQC(reads)
    FASTP(reads)
    POST_FASTQC(FASTP.out.R1_trimmed, FASTP.out.R2_trimmed)
    MULTIQC()
}
