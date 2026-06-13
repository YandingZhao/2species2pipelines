#!/usr/bin/env python3
"""Run OrthoFinder and extract one-to-one orthologs as a gene map TSV.

OrthoFinder must be installed and on PATH (conda install -c bioconda orthofinder).

Output: gene_map.tsv (two columns, no header):
  col1 = species_a gene symbol
  col2 = species_b gene symbol

This file is consumed by run_ortholog_convert.R --strategy orthofinder.

Usage
-----
  python3 run_orthofinder.py \\
    --fasta_a  dog_proteins.fa \\
    --fasta_b  human_proteins.fa \\
    --species_a dog \\
    --species_b human \\
    --threads   16 \\
    --outdir    tests/data/ortholog_maps/orthofinder
"""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="OrthoFinder one-to-one ortholog map")
    p.add_argument("--fasta_a",   required=True)
    p.add_argument("--fasta_b",   required=True)
    p.add_argument("--species_a", required=True)
    p.add_argument("--species_b", required=True)
    p.add_argument("--threads",   type=int, default=8)
    p.add_argument("--outdir",    required=True)
    return p.parse_args()


def _check_orthofinder():
    if shutil.which("orthofinder") is None:
        raise RuntimeError(
            "OrthoFinder not found on PATH. "
            "Install with: conda install -c bioconda orthofinder"
        )


def _run_orthofinder(fasta_dir: str, threads: int, outdir: str) -> Path:
    cmd = [
        "orthofinder",
        "-f",  fasta_dir,
        "-t",  str(threads),
        "-o",  outdir,
        "-og",            # stop after orthogroup inference (faster, skip gene trees)
    ]
    subprocess.run(cmd, check=True)
    # OrthoFinder creates a timestamped results dir inside outdir
    results_dirs = sorted(Path(outdir).glob("Results_*"))
    if not results_dirs:
        raise RuntimeError("OrthoFinder produced no Results_* directory")
    return results_dirs[-1]


def _parse_orthogroups(results_dir: Path, species_a: str, species_b: str) -> list:
    """Extract 1-to-1 orthologs from OrthoFinder's Orthogroups.tsv."""
    og_file = results_dir / "Orthogroups" / "Orthogroups.tsv"
    if not og_file.exists():
        raise FileNotFoundError(f"Expected {og_file}")

    with open(og_file) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Columns: Orthogroup | species1_name | species2_name | ...
        # Match by species name (filename stem, stripped of .fa/.fasta)
        col_a = col_b = None
        for i, h in enumerate(header[1:], 1):
            stem = Path(h).stem if "." in h else h
            if species_a.lower() in stem.lower():
                col_a = i
            if species_b.lower() in stem.lower():
                col_b = i
        if col_a is None or col_b is None:
            raise RuntimeError(
                f"Could not find species columns in {og_file}.\n"
                f"Header: {header}\nSpecies: {species_a}, {species_b}"
            )

        pairs = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(col_a, col_b):
                continue
            genes_a = [g.strip() for g in parts[col_a].split(",") if g.strip()]
            genes_b = [g.strip() for g in parts[col_b].split(",") if g.strip()]
            # Strict 1-to-1 only
            if len(genes_a) == 1 and len(genes_b) == 1:
                pairs.append((genes_a[0], genes_b[0]))

    return pairs


def main():
    args = parse_args()
    _check_orthofinder()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        # OrthoFinder reads all FASTAs from a single directory
        fasta_dir = os.path.join(tmpdir, "fastas")
        os.makedirs(fasta_dir)
        ext_a = Path(args.fasta_a).suffix
        ext_b = Path(args.fasta_b).suffix
        shutil.copy(args.fasta_a, os.path.join(fasta_dir, f"{args.species_a}{ext_a}"))
        shutil.copy(args.fasta_b, os.path.join(fasta_dir, f"{args.species_b}{ext_b}"))

        of_outdir = os.path.join(tmpdir, "orthofinder_out")
        print("Running OrthoFinder...")
        results_dir = _run_orthofinder(fasta_dir, args.threads, of_outdir)
        print(f"OrthoFinder results: {results_dir}")

        pairs = _parse_orthogroups(results_dir, args.species_a, args.species_b)

    out_map = outdir / "gene_map.tsv"
    with open(out_map, "w") as fh:
        for g_a, g_b in pairs:
            fh.write(f"{g_a}\t{g_b}\n")

    out_report = outdir / "orthofinder_report.txt"
    with open(out_report, "w") as fh:
        fh.write(f"species_a: {args.species_a}\n")
        fh.write(f"species_b: {args.species_b}\n")
        fh.write(f"one_to_one_pairs: {len(pairs)}\n")
        fh.write("status: ok\n")

    print(f"One-to-one ortholog pairs: {len(pairs)}")
    print(f"Map written to: {out_map}")


if __name__ == "__main__":
    main()
