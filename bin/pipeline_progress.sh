#!/usr/bin/env bash
# Percentage of genomes that have finished each pipeline step.
#
# Usage: bash bin/pipeline_progress.sh <samplesheet.tsv> [results_dir]
# Writes:  <results_dir>/progress.txt             summary table (done / failed / pending %)
#          <results_dir>/progress_per_sample.tsv  one row per genome: done | FAILED | -
# A step counts as done when its final output file exists in the genome's folder.
# Failures come from Nextflow's trace file (last attempt failed, errorStrategy 'ignore').
set -u

sheet=${1:?usage: pipeline_progress.sh <samplesheet.tsv> [results_dir]}
outdir=${2:-results}
trace=$outdir/pipeline_info/trace.txt
steps=(fastp spades quast checkm2 bakta gtdbtk)

expected() {  # <step> <sample> -> file that marks the step as finished
    case $1 in
        fastp)   echo "$outdir/$2/fastp/$2.fastp.json" ;;
        spades)  echo "$outdir/$2/spades/$2.fasta" ;;
        quast)   echo "$outdir/$2/quast/report.tsv" ;;
        checkm2) echo "$outdir/$2/checkm2/quality_report.tsv" ;;
        bakta)   echo "$outdir/$2/bakta/$2.gff3" ;;
        gtdbtk)  echo "$outdir/$2/gtdbtk/$2.gtdbtk.summary.tsv" ;;
    esac
}

mapfile -t samples < <(tr -d '\r' < "$sheet" | awk -F'\t' 'NF >= 3 && $1 != "" && $1 !~ /^#/ && tolower($1) != "sample_name" { print $1 }')
total=${#samples[@]}
mkdir -p "$outdir"

# Failed tasks: "<step>|<sample>". GTDB-Tk runs per batch, so its staged genomes give the samples.
# n_failed is kept separately: with set -u, bash treats an empty associative array as unset
declare -A failed
n_failed=0
if [[ -s $trace ]]; then
    while IFS=$'\t' read -r proc tag wdir; do
        step=${proc,,}
        if [[ $step == gtdbtk ]]; then
            for fa in "$wdir"/genomes/*.fasta; do
                [[ -e $fa ]] && failed["gtdbtk|$(basename "$fa" .fasta)"]=1 && n_failed=$(( n_failed + 1 ))
            done
        else
            failed["$step|$tag"]=1
            n_failed=$(( n_failed + 1 ))
        fi
    done < <(awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        { k = $c["process"] "\t" $c["tag"]
          if ($c["status"] == "COMPLETED" || $c["status"] == "CACHED") ok[k] = 1
          else if ($c["status"] == "FAILED" || $c["status"] == "ABORTED") { bad[k] = 1; wd[k] = $c["workdir"] } }
        END { for (k in bad) if (!(k in ok)) print k "\t" wd[k] }' "$trace")
fi

declare -A n_done n_fail
for st in "${steps[@]}"; do n_done[$st]=0; n_fail[$st]=0; done
all_done=0

{
    printf 'sample'; printf '\t%s' "${steps[@]}"; printf '\n'
    for s in "${samples[@]}"; do
        printf '%s' "$s"
        complete=1
        for st in "${steps[@]}"; do
            if [[ -s $(expected "$st" "$s") ]]; then
                status=done; n_done[$st]=$(( n_done[$st] + 1 ))
            elif [[ -n ${failed["$st|$s"]:-} ]]; then
                status=FAILED; n_fail[$st]=$(( n_fail[$st] + 1 )); complete=0
            else
                status=-; complete=0
            fi
            printf '\t%s' "$status"
        done
        printf '\n'
        (( complete )) && all_done=$(( all_done + 1 ))
    done
} > "$outdir/progress_per_sample.tsv"

pct() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%5.1f%%", b ? 100 * a / b : 0 }'; }
bar() { awk -v a="$1" -v b="$2" 'BEGIN { n = b ? int(20 * a / b) : 0; s = ""; for (i = 0; i < 20; i++) s = s (i < n ? "#" : "."); print s }'; }

{
    echo "Genome pipeline progress -- updated $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Samplesheet: $sheet"
    echo "Genomes:     $total"
    echo
    printf '%-10s %7s %7s %7s %8s  %s\n' step done failed pending percent progress
    for st in "${steps[@]}"; do
        d=${n_done[$st]}; f=${n_fail[$st]}
        printf '%-10s %7d %7d %7d %8s  [%s]\n' "$st" "$d" "$f" $(( total - d - f )) "$(pct "$d" "$total")" "$(bar "$d" "$total")"
    done
    printf '%-10s %7d %7s %7s %8s  [%s]\n' "ALL STEPS" "$all_done" "" "" "$(pct "$all_done" "$total")" "$(bar "$all_done" "$total")"
    echo
    echo "Per-genome status: $outdir/progress_per_sample.tsv"
    if (( n_failed > 0 )); then
        echo
        echo "Failed (see the task's .command.err in the work dir listed in $trace):"
        for k in "${!failed[@]}"; do echo "  ${k%%|*}  ${k#*|}"; done | sort
    fi
} > "$outdir/progress.txt"

cat "$outdir/progress.txt"
