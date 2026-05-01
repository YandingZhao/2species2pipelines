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
- `scripts/run_harmony_module.R` as the Harmony runner inspired by benchmark scripts
- `scripts/run_seurat4_module.R` as the Seurat4 CCA runner inspired by benchmark scripts
- `scripts/run_fastmnn_module.R` as the fastMNN runner inspired by benchmark scripts
- `scripts/run_bbknn_module.py` as the BBKNN runner inspired by benchmark scripts
- `scripts/generate_bbknn_test_data.py` to generate synthetic `.h5ad` test data
- `docker/Dockerfile` as the integration runtime image
- `docker/Dockerfile.bbknn` as the BBKNN runtime image
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
docker run --rm -v "$PWD":/work -w /work local/bbknn-module:dev python scripts/generate_bbknn_test_data.py
```

```bash
nextflow run . -profile test,docker -stub-run \
	--harmony_container local/harmony-module:dev \
	--seurat4_container local/harmony-module:dev \
	--fastmnn_container local/harmony-module:dev \
	--bbknn_container local/bbknn-module:dev
```

Run integration modules without stub:

```bash
nextflow run . -profile docker,test \
	--harmony_container local/harmony-module:dev \
	--seurat4_container local/harmony-module:dev \
	--fastmnn_container local/harmony-module:dev \
	--bbknn_container local/bbknn-module:dev
```

Outputs are written to `results/` (or `tests/results/` with the test profile).

## CI

On push to `main` and pull requests, GitHub Actions will:
- build `docker/Dockerfile`
- run `nextflow run . -profile test,docker -stub-run`