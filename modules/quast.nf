process QUAST {
    tag "$sample"
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    path "quast", emit: dir

    script:
    """
    quast.py ${fasta} \\
        -o quast \\
        -t ${task.cpus} \\
        --labels ${sample} \\
        --min-contig ${params.min_contig_len}
    """

    stub:
    """
    mkdir -p quast
    printf 'Assembly\\tN50\\n${sample}\\t1000\\n' > quast/transposed_report.tsv
    cp quast/transposed_report.tsv quast/report.tsv
    """
}
