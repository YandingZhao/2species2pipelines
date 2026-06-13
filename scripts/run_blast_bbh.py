#!/usr/bin/env python3
"""Compute BLAST bidirectional best hits (BBH) between two species' proteomes.

Outputs a two-column TSV (no header):
  col1 = species_a gene symbol
  col2 = species_b gene symbol

This file is consumed by run_ortholog_convert.R --strategy blast_bbh_XX.

Prerequisites
-------------
  BLAST+ (blastp) must be on PATH.
  Input FASTA files: one protein sequence per gene, FASTA header = gene symbol.

Typical usage
-------------
  # Download protein FASTAs from Ensembl BioMart first, then:
  python3 run_blast_bbh.py \\
    --fasta_a  dog_proteins.fa \\
    --fasta_b  human_proteins.fa \\
    --identity 80 \\
    --threads  8 \\
    --outdir   tests/data/ortholog_maps/blast_bbh_80
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="BLAST bidirectional best hits")
    p.add_argument("--fasta_a",   required=True, help="Species A protein FASTA")
    p.add_argument("--fasta_b",   required=True, help="Species B protein FASTA")
    p.add_argument("--identity",  type=float, default=80.0,
                   help="Minimum %% amino acid identity threshold (default: 80)")
    p.add_argument("--coverage",  type=float, default=50.0,
                   help="Minimum query coverage %% (default: 50)")
    p.add_argument("--evalue",    type=float, default=1e-5,
                   help="Maximum e-value (default: 1e-5)")
    p.add_argument("--threads",   type=int,   default=4)
    p.add_argument("--outdir",    required=True,
                   help="Output directory (will be created if needed)")
    return p.parse_args()


def make_blast_db(fasta: str, db_path: str) -> None:
    cmd = ["makeblastdb", "-in", fasta, "-dbtype", "prot", "-out", db_path]
    subprocess.run(cmd, check=True, capture_output=True)


def run_blastp(query: str, db: str, out: str, threads: int, evalue: float) -> None:
    cmd = [
        "blastp",
        "-query",   query,
        "-db",      db,
        "-out",     out,
        "-outfmt",  "6 qseqid sseqid pident length qlen slen evalue bitscore",
        "-evalue",  str(evalue),
        "-num_threads", str(threads),
        "-max_target_seqs", "1",   # only top hit per query
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def parse_blast(blast_file: str, min_identity: float, min_coverage: float) -> dict:
    """Return dict: query_gene -> (subject_gene, pident)  (only top hit per query)."""
    hits = {}
    with open(blast_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            qid, sid, pident, length, qlen, slen, evalue, bitscore = parts
            pident = float(pident)
            length = int(length)
            qlen   = int(qlen)
            cov    = 100.0 * length / qlen if qlen > 0 else 0.0

            if pident < min_identity or cov < min_coverage:
                continue
            # keep best (first, since blastp output is sorted by bitscore)
            if qid not in hits:
                hits[qid] = sid
    return hits


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        db_a = os.path.join(tmpdir, "db_a")
        db_b = os.path.join(tmpdir, "db_b")
        blast_atob = os.path.join(tmpdir, "a_vs_b.txt")
        blast_btoa = os.path.join(tmpdir, "b_vs_a.txt")

        print("Building BLAST databases...")
        make_blast_db(args.fasta_a, db_a)
        make_blast_db(args.fasta_b, db_b)

        print("Running BLAST A→B...")
        run_blastp(args.fasta_a, db_b, blast_atob, args.threads, args.evalue)

        print("Running BLAST B→A...")
        run_blastp(args.fasta_b, db_a, blast_btoa, args.threads, args.evalue)

    print("Parsing hits and finding bidirectional best hits...")
    hits_atob = parse_blast(blast_atob, args.identity, args.coverage)
    hits_btoa = parse_blast(blast_btoa, args.identity, args.coverage)

    # BBH: gene_a → gene_b AND gene_b → gene_a (reciprocal)
    bbh_pairs = []
    for gene_a, gene_b in hits_atob.items():
        if hits_btoa.get(gene_b) == gene_a:
            bbh_pairs.append((gene_a, gene_b))

    out_map = outdir / "gene_map.tsv"
    with open(out_map, "w") as fh:
        for gene_a, gene_b in bbh_pairs:
            fh.write(f"{gene_a}\t{gene_b}\n")

    out_report = outdir / "blast_bbh_report.txt"
    with open(out_report, "w") as fh:
        fh.write(f"identity_threshold: {args.identity}\n")
        fh.write(f"coverage_threshold: {args.coverage}\n")
        fh.write(f"evalue_threshold: {args.evalue}\n")
        fh.write(f"a_to_b_hits: {len(hits_atob)}\n")
        fh.write(f"b_to_a_hits: {len(hits_btoa)}\n")
        fh.write(f"bbh_pairs: {len(bbh_pairs)}\n")
        fh.write("status: ok\n")

    print(f"BBH pairs found: {len(bbh_pairs)}")
    print(f"Map written to: {out_map}")


if __name__ == "__main__":
    main()
