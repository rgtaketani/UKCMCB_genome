#!/usr/bin/env nextflow

include { FASTP   } from './modules/fastp'
include { SPADES  } from './modules/spades'
include { QUAST   } from './modules/quast'
include { CHECKM2 } from './modules/checkm2'
include { BAKTA   } from './modules/bakta'
include { GTDBTK  } from './modules/gtdbtk'

// The three databases must be given and exist; checked before anything is submitted.
def checkDatabases() {
    def dbs = [bakta_db: params.bakta_db, checkm2_db: params.checkm2_db, gtdbtk_db: params.gtdbtk_db]
    def errors = []
    dbs.each { name, value ->
        if (!value) {
            errors << "--${name} is not set"
        }
        else if (!file(value).exists()) {
            errors << "--${name} not found: ${value}"
        }
    }
    if (errors) {
        error("Reference databases missing (see README, 'Databases'):\n  " + errors.join('\n  '))
    }
}

def readPath(p, sheet_dir) {
    return p.startsWith('/') || p.contains('://') ? file(p) : sheet_dir.resolve(p)
}

// Samplesheet: sample_name <tab> read1 <tab> read2. A header line starting with
// "sample_name" and lines starting with '#' are skipped. Relative read paths are taken
// relative to the samplesheet's folder. All problems are reported at once.
def parseSamplesheet(path) {
    def rows = []
    def errors = []
    def seen = [:]
    def sheet_dir = file(path).parent
    file(path).readLines().eachWithIndex { raw, i ->
        def line = raw.replaceAll('[\\uFEFF\\r]', '')
        def cols = line.split('\t').collect { c -> c.trim() }
        def skip = !line.trim() || line.startsWith('#') || cols[0].toLowerCase() == 'sample_name'
        if (!skip && cols.size() < 3) {
            errors << "line ${i + 1}: expected 3 tab-separated columns, found ${cols.size()}"
        }
        else if (!skip) {
            def name = cols[0]
            def r1 = readPath(cols[1], sheet_dir)
            def r2 = readPath(cols[2], sheet_dir)
            if (!(name ==~ /[A-Za-z0-9._-]+/)) {
                errors << "line ${i + 1}: sample name '${name}' may only contain letters, digits, '.', '_' and '-'"
            }
            if (seen.containsKey(name)) {
                errors << "line ${i + 1}: duplicate sample name '${name}'"
            }
            seen[name] = true
            [r1, r2].each { r ->
                if (!r.exists()) {
                    errors << "line ${i + 1}: read file not found: ${r}"
                }
            }
            rows << [name, r1, r2]
        }
    }
    if (errors) error("Samplesheet ${path} has ${errors.size()} problem(s):\n  " + errors.join('\n  '))
    if (!rows) error("Samplesheet ${path} contains no samples")
    return rows
}

workflow {
    main:
    if (!params.samplesheet) error("Please provide --samplesheet <file.tsv>")
    checkDatabases()

    def samples = parseSamplesheet(params.samplesheet)

    // GTDB-Tk batches are fixed by samplesheet order so -resume reuses finished batches
    def bs = params.gtdbtk_batch_size as int
    def batch_of = [:]
    samples.eachWithIndex { s, i -> batch_of[s[0]] = String.format('%03d', i.intdiv(bs) + 1) }
    def batch_size = batch_of.values().countBy { b -> b }

    log.info """
    ==============================================
     Isolate genome pipeline
     samplesheet : ${params.samplesheet}
     genomes     : ${samples.size()}
     GTDB batches: ${batch_size.size()} (up to ${bs} genomes each)
     outdir      : ${params.outdir}
    ==============================================
    """.stripIndent()

    FASTP(channel.fromList(samples))
    SPADES(FASTP.out.reads)

    QUAST(SPADES.out.assembly)
    CHECKM2(SPADES.out.assembly)
    BAKTA(SPADES.out.assembly)

    // A batch starts as soon as all its genomes are assembled; batches with failed
    // assemblies are released (incomplete) when assembly finishes for every genome
    def gtdb_in = SPADES.out.assembly
        .map { sample, fasta -> tuple(groupKey(batch_of[sample], batch_size[batch_of[sample]]), fasta) }
        .groupTuple(remainder: true, sort: true)
        .map { key, fastas -> tuple(key.toString(), fastas) }
    GTDBTK(gtdb_in)

    // Combined tables across all genomes
    QUAST.out.dir.map { d -> d.resolve('transposed_report.tsv') }
        .collectFile(name: 'quast_summary.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/summary")
    CHECKM2.out.dir.map { d -> d.resolve('quality_report.tsv') }
        .collectFile(name: 'checkm2_summary.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/summary")
    GTDBTK.out.summaries.flatten()
        .collectFile(name: 'gtdbtk_summary.tsv', keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/summary")

    onComplete:
    log.info "Pipeline finished: ${workflow.success ? 'OK' : 'with errors'} -- see ${params.outdir}/progress.txt for per-step status"
}
