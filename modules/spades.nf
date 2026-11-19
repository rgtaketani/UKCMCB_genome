process SPADES {
    tag "$sample"
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("spades/${sample}.fasta"), emit: assembly
    path "spades/${sample}.spades.*"

    script:
    """
    spades.py ${params.spades_mode} \\
        -1 ${r1} -2 ${r2} \\
        -o spades_run \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()}

    mkdir -p spades
    cp spades_run/contigs.fasta                       spades/${sample}.spades.contigs.fasta
    cp spades_run/spades.log                          spades/${sample}.spades.log
    cp spades_run/scaffolds.fasta                     spades/${sample}.spades.scaffolds.fasta      2>/dev/null || true
    cp spades_run/assembly_graph_with_scaffolds.gfa   spades/${sample}.spades.assembly_graph.gfa   2>/dev/null || true

    # Final assembly used downstream: contigs >= min_contig_len (length is in the SPAdes header)
    awk -v min=${params.min_contig_len} '/^>/ { split(\$0, a, "_"); keep = (a[4] + 0 >= min) } keep' \\
        spades_run/contigs.fasta > spades/${sample}.fasta

    if [ ! -s spades/${sample}.fasta ]; then
        echo "ERROR: no contigs >= ${params.min_contig_len} bp for ${sample}" >&2
        exit 1
    fi

    # Intermediate k-mer folders and corrected reads are large; drop them from the work dir
    rm -rf spades_run/K* spades_run/corrected spades_run/tmp spades_run/misc
    """

    stub:
    """
    mkdir -p spades
    printf '>NODE_1_length_1000_cov_10\\nACGT\\n' > spades/${sample}.fasta
    cp spades/${sample}.fasta spades/${sample}.spades.contigs.fasta
    touch spades/${sample}.spades.log
    """
}
