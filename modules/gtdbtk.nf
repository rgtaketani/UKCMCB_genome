// GTDB-Tk loads the reference data once per job, so genomes are classified in batches;
// the result for each genome is then written into its own folder.
process GTDBTK {
    tag "batch_${batch_id}"
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    tuple val(batch_id), path(fastas, stageAs: 'genomes/*')

    output:
    path "gtdbtk_batches/batch_${batch_id}", emit: batch_dir
    path "*/gtdbtk/*.gtdbtk.summary.tsv", emit: summaries

    script:
    def full_tree = params.gtdbtk_full_tree ? '--full_tree' : ''
    """
    export GTDBTK_DATA_PATH=${params.gtdbtk_db}
    od=gtdbtk_batches/batch_${batch_id}
    rc=0
    gtdbtk classify_wf \\
        --genome_dir genomes \\
        --out_dir \$od \\
        -x fasta \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${params.gtdbtk_pplacer_cpus} \\
        ${full_tree} || rc=\$?

    # GTDB-Tk can exit non-zero when every genome is filtered (too few markers) or
    # during temp-dir cleanup; only treat it as a failure if it produced nothing
    if [ \$rc -ne 0 ] && ! ls \$od/gtdbtk.*.summary.tsv >/dev/null 2>&1 && [ ! -s \$od/gtdbtk.filtered.tsv ]; then
        exit \$rc
    fi

    hdr=\$(cat \$od/gtdbtk.bac120.summary.tsv \$od/gtdbtk.ar53.summary.tsv 2>/dev/null | head -n 1)
    [ -n "\$hdr" ] || hdr=\$(printf 'user_genome\\tclassification\\tnote')
    ncol=\$(printf '%s\\n' "\$hdr" | awk -F'\\t' '{ print NF }')

    for fa in genomes/*.fasta; do
        s=\$(basename "\$fa" .fasta)
        mkdir -p "\$s/gtdbtk"
        out="\$s/gtdbtk/\$s.gtdbtk.summary.tsv"
        row=\$(awk -F'\\t' -v s="\$s" '\$1 == s' \$od/gtdbtk.*.summary.tsv 2>/dev/null | head -n 1)
        if [ -z "\$row" ]; then
            # Not classified: keep the same columns, put the reason in the last one
            why=\$(awk -F'\\t' -v s="\$s" '\$1 == s { print \$2 }' \$od/gtdbtk.filtered.tsv \$od/gtdbtk.failed_genomes.tsv 2>/dev/null | head -n 1)
            row=\$(awk -v s="\$s" -v n="\$ncol" -v why="\${why:-not classified by GTDB-Tk}" \\
                'BEGIN { r = s "\\tUnclassified"; for (i = 3; i < n; i++) r = r "\\tN/A"; print r "\\t" why }')
        fi
        printf '%s\\n%s\\n' "\$hdr" "\$row" > "\$out"
    done
    """

    stub:
    """
    mkdir -p gtdbtk_batches/batch_${batch_id}
    for fa in genomes/*.fasta; do
        s=\$(basename "\$fa" .fasta)
        mkdir -p "\$s/gtdbtk"
        printf 'user_genome\\tclassification\\n%s\\td__Bacteria\\n' "\$s" > "\$s/gtdbtk/\$s.gtdbtk.summary.tsv"
    done
    """
}
