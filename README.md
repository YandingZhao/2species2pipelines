# species to pipelines 2

This project is an attempt to recreate the pipeline located in the benchmark_project_script. While recreating, we do not want to use the same code, but instead want to create a more modular and reusable codebase. The goal is to create a pipeline that can be easily for integrating canine and human data.

## Tech stack
- Use nextflow as the main pipeline engine
- Use nf-test for pipeline tests
- Use docker for containerization

## nf-core base pipeline scaffold

This repository now includes a minimal nf-core style Nextflow scaffold.

Key files:
- `main.nf` as the workflow entrypoint
- `workflows/nfcore_base.nf` as the base workflow definition
- `modules/local/make_run_metadata.nf` as a local example process
- `modules/local/harmony_integration.nf` as the first integration module
- `modules/local/seurat4_integration.nf` as the second integration module
- `modules/local/fastmnn_integration.nf` as the third integration module
- `modules/local/bbknn_integration.nf` as the first Python integration module
- `modules/local/scanorama_integration.nf` as the second Python integration module
- `modules/local/ortholog_convert_pair.nf` as the first preprocessing step for cross-species gene harmonization
- `modules/local/seurat_to_anndata_pair.nf` for automatic Seurat `.rds` to `.h5ad` conversion for BBKNN inputs
- `scripts/run_harmony_module.R` as the Harmony runner inspired by benchmark scripts
- `scripts/run_seurat4_module.R` as the Seurat4 CCA runner inspired by benchmark scripts
- `scripts/run_fastmnn_module.R` as the fastMNN runner inspired by benchmark scripts
- `scripts/run_bbknn_module.py` as the BBKNN runner inspired by benchmark scripts
- `scripts/run_scanorama_module.py` as the Scanorama runner inspired by benchmark scripts
- `scripts/run_ortholog_convert_pair.R` for species_a -> species_b ortholog conversion on Seurat `.rds` inputs
- `scripts/run_seurat_to_anndata_pair.R` for converting Seurat pair inputs to `.h5ad`
- `docker/Dockerfile` as the integration runtime image
- `docker/Dockerfile.bbknn` as the shared Python runtime image for BBKNN and Scanorama
- `.github/workflows/docker-and-nextflow.yml` as CI build and smoke test
- `nextflow.config` with `standard`, `test`, and `docker` profiles
- `conf/base.config` and `conf/test.config` for profile-specific settings
- `assets/samplesheet.csv` as the default input
- `params.schema.json` for parameter schema metadata

Run locally:

```bash
nextflow run . -profile test -stub-run
```

Run with docker (build image first):

```bash
docker build -t local/harmony-module:dev -f docker/Dockerfile .
```

```bash
docker build -t local/bbknn-module:dev -f docker/Dockerfile.bbknn .
```

```bash
nextflow run . -profile test,docker -stub-run \
	--harmony_container local/harmony-module:dev \
	--seurat4_container local/harmony-module:dev \
	--fastmnn_container local/harmony-module:dev \
	--bbknn_container local/bbknn-module:dev \
	--scanorama_container local/bbknn-module:dev
```

Run integration modules without stub:

```bash
nextflow run . -profile docker,test \
	--harmony_container local/harmony-module:dev \
	--seurat4_container local/harmony-module:dev \
	--fastmnn_container local/harmony-module:dev \
	--bbknn_container local/bbknn-module:dev \
	--scanorama_container local/bbknn-module:dev
```

Outputs are written to `results/` (or `tests/results/` with the test profile).

For `.rds` inputs, ortholog conversion is the first step for all integration modules.
Gene symbols from `species_a` are converted to `species_b` with `orthogene` (`non121_strategy=drop_both_species`, `method=gprofiler`), then both inputs are restricted to shared genes.

Use species names supported by `orthogene` in the samplesheet (for example `dog` instead of `canine`).

For BBKNN, if `source_a` and `source_b` are Seurat `.rds` files in the input samplesheet, they are converted automatically to `.h5ad` before BBKNN runs.

For Scanorama, the same `.rds` to `.h5ad` conversion path is used automatically before integration.

## CI

On push to `main` and pull requests, GitHub Actions will:
- build `docker/Dockerfile`
- run `nextflow run . -profile test,docker -stub-run`