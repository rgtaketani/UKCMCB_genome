process CHECKM2 {
    tag "$sample"
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    path "checkm2", emit: dir

    script:
    """
    checkm2 predict \\
        --input ${fasta} \\
        --output-directory checkm2 \\
        --threads ${task.cpus} \\
        --database_path ${params.checkm2_db}
    """

    stub:
    """
    mkdir -p checkm2
    printf 'Name\\tCompleteness\\tContamination\\n${sample}\\t99.0\\t0.5\\n' > checkm2/quality_report.tsv
    """
}
