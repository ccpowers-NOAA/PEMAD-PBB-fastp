process pre_fastqc {
    executor params.pre_fastqc.executor

    publishDir "$launchDir/${params.publishDir}", overwrite: true, pattern: "${id}/annotate/NF_work_dir_annotate.txt", mode: 'copy'

    errorStrategy 'finish'

    tag "${id}"

    input:
        tuple val(id), val(path), path(assembly), path(coverage), val(opts)

    output:
    tuple val(id), val(path),
        path("${id}/annotate/${id}_annotations_*.csv"),
        path("${id}/annotate/${id}_assembly_*.fasta"),
        path("${id}/annotate/${id}_coverageStats_*.csv"),
        path("${id}/annotate/NF_work_dir_annotate.txt")                 // Nextflow working directory, for troubleshooting


    shell:
    '''
    mkdir -p log
    fastqc 00_reads/!{sample}
    filename=!{sample}
    fastqc 00_reads/${{filename/R1.fastq.gz/R2.fastq.gz}}
    touch !{output}
    '''
}

workflow pre_fastqc_wf {
    channel.fromQuery(params.sqlRead, db: 'sqlite')
        .map{ it ->
            tuple(
                it[0],                                          // ID
                it[1],                                          // path
                file(                                           // Assembly
                    params.publishDir + '/' +
                    it[0] + '/assemble/' + it[2] + '/' +
                    it[0] + '_assembly_' + it[1] + '.fasta'
                ),
                file(                                           // Coverage
                    params.publishDir + '/' +
                    it[0] + '/assemble/' + it[2] + '/' +
                    it[0] + '_assembly_' + it[1] + '_coverageStats.csv'
                ),
                [
                    cpus:  it[3],                                      // cpus
                    memory: it[4],                                     // memory
                    ref_db: it[5],                                     // mitos_ref_db
                    ref_dir: it[6],                                    // mitos_ref_dir
                    mitos: it[7],                                      // mitos_opts
                    trnaScan: it[8],                                    // trnaScan_opts
                    start_gene: it[9]                                  // starting gene for rotation
                ]

            )
        }
        .set { annotate_in }

    pre_fastqc(annotate_in).set { annotate_out }

    emit:
           ch = annotate_out[0]

}
