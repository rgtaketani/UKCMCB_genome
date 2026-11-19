#!/usr/bin/env python3
"""Simulate Illumina-like paired-end reads (gzipped FASTQ) from a reference FASTA.

Used to build the small diagnostic dataset in test/reads/. No dependencies.

Usage: simulate_reads.py <ref.fasta> <out_prefix> <n_pairs> [seed]
Writes <out_prefix>_R1.fastq.gz and <out_prefix>_R2.fastq.gz
"""
import gzip
import random
import sys

READ_LEN = 150
INSERT_MEAN, INSERT_SD = 350, 40
SUB_RATE = 0.002
# TruSeq adapter, appended when the insert is shorter than the read (tests fastp trimming)
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def mutate(read, rng):
    bases = list(read)
    for i in range(len(bases)):
        if rng.random() < SUB_RATE:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def quals(n, rng):
    # High quality that decays towards the 3' end
    return "".join(chr(33 + max(2, min(40, int(rng.gauss(38 - 10 * i / n, 3))))) for i in range(n))


def main():
    ref, prefix, n_pairs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    genome = read_fasta(ref)
    glen = len(genome)
    with gzip.open(f"{prefix}_R1.fastq.gz", "wt", compresslevel=6) as r1, \
         gzip.open(f"{prefix}_R2.fastq.gz", "wt", compresslevel=6) as r2:
        for i in range(n_pairs):
            ins = max(80, int(rng.gauss(INSERT_MEAN, INSERT_SD)))
            start = rng.randrange(0, glen - ins)
            frag = genome[start:start + ins]
            if rng.random() < 0.5:
                frag = frag.translate(COMP)[::-1]
            s1 = (frag + ADAPTER_R1 + "A" * READ_LEN)[:READ_LEN]
            s2 = (frag.translate(COMP)[::-1] + ADAPTER_R2 + "A" * READ_LEN)[:READ_LEN]
            name = f"@sim_{i + 1}"
            r1.write(f"{name} 1:N:0:1\n{mutate(s1, rng)}\n+\n{quals(READ_LEN, rng)}\n")
            r2.write(f"{name} 2:N:0:1\n{mutate(s2, rng)}\n+\n{quals(READ_LEN, rng)}\n")


if __name__ == "__main__":
    main()
