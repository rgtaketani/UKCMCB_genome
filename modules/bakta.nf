process BAKTA {
    tag "$sample"
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    path "bakta", emit: dir

    script:
    """
    bakta --db ${params.bakta_db} \\
        --output bakta \\
        --prefix ${sample} \\
        --threads ${task.cpus} \\
        --min-contig-length ${params.min_contig_len} \\
        --tmp-dir \$PWD \\
        ${fasta}
    """

    stub:
    """
    mkdir -p bakta
    echo '##gff-version 3' > bakta/${sample}.gff3
    """
}
