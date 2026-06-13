# Ortholog Maps

Pre-computed gene mapping tables for BLAST BBH and OrthoFinder strategies.

## Format

Each `gene_map.tsv` is a two-column TSV (no header):
```
col1 = species_a gene symbol
col2 = species_b gene symbol
```

## How to Generate

### BLAST BBH maps

```bash
# Download protein FASTAs from Ensembl BioMart (one entry per gene symbol)
# Then run:
python3 scripts/run_blast_bbh.py \
  --fasta_a  dog_proteins.fa \
  --fasta_b  human_proteins.fa \
  --identity 60 \
  --threads  8 \
  --outdir   tests/data/ortholog_maps/blast_bbh_60

python3 scripts/run_blast_bbh.py \
  --fasta_a  dog_proteins.fa \
  --fasta_b  human_proteins.fa \
  --identity 80 \
  --outdir   tests/data/ortholog_maps/blast_bbh_80

python3 scripts/run_blast_bbh.py \
  --fasta_a  dog_proteins.fa \
  --fasta_b  human_proteins.fa \
  --identity 90 \
  --outdir   tests/data/ortholog_maps/blast_bbh_90
```

### OrthoFinder map

```bash
python3 scripts/run_orthofinder.py \
  --fasta_a   dog_proteins.fa \
  --fasta_b   human_proteins.fa \
  --species_a dog \
  --species_b human \
  --threads   16 \
  --outdir    tests/data/ortholog_maps/orthofinder
```

Commit the resulting `gene_map.tsv` files to enable those strategies in CI.
