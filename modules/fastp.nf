process FASTP {
    tag "$sample"
    publishDir path: { "${params.outdir}/${sample}" }, mode: 'copy',
        saveAs: { fn -> !params.save_trimmed_reads && fn.endsWith('.fastq.gz') ? null : fn }

    input:
    tuple val(sample), path(r1, stageAs: 'raw/R1.fastq.gz'), path(r2, stageAs: 'raw/R2.fastq.gz')

    output:
    tuple val(sample), path("fastp/${sample}_R1.trimmed.fastq.gz"), path("fastp/${sample}_R2.trimmed.fastq.gz"), emit: reads
    path "fastp/${sample}.fastp.{json,html,log}"

    script:
    """
    mkdir -p fastp
    fastp \\
        -i ${r1} -I ${r2} \\
        -o fastp/${sample}_R1.trimmed.fastq.gz \\
        -O fastp/${sample}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --json fastp/${sample}.fastp.json \\
        --html fastp/${sample}.fastp.html \\
        2> fastp/${sample}.fastp.log
    """

    stub:
    """
    mkdir -p fastp
    touch fastp/${sample}_R1.trimmed.fastq.gz fastp/${sample}_R2.trimmed.fastq.gz
    echo '{}' > fastp/${sample}.fastp.json
    touch fastp/${sample}.fastp.html fastp/${sample}.fastp.log
    """
}
