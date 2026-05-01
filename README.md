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
- `scripts/run_harmony_module.R` as the Harmony runner inspired by benchmark scripts
- `docker/harmony/Dockerfile` as the Harmony runtime image
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
docker build -t local/harmony-module:dev -f docker/harmony/Dockerfile .
```

```bash
nextflow run . -profile test,docker -stub-run --harmony_container local/harmony-module:dev
```

Run Harmony module without stub:

```bash
nextflow run . -profile docker,test --harmony_container local/harmony-module:dev
```

Outputs are written to `results/` (or `tests/results/` with the test profile).

## CI

On push to `main` and pull requests, GitHub Actions will:
- build `docker/harmony/Dockerfile`
- run `nextflow run . -profile test,docker -stub-run`