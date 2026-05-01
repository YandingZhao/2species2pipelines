# Project Guidelines

## Architecture
- The current codebase is script-based and lives under `benchmark_project_script/`, it's not runnable.
- The root project goal is a future modular Nextflow rewrite (see `README.md`), but there are currently no Nextflow pipeline files in this repository.
- Shared method/evaluation logic is in `benchmark_project_script/core_script/`; method runners and evaluation entrypoints are in `benchmark_project_script/`.
- Keep changes focused on existing Python/R benchmark scripts unless the task explicitly asks for Nextflow scaffolding.
- Minified test data as located in /test_data/ should be used for any testing or development work to avoid large file handling issues.

## Build and Test
- Treat benchmark_project_script/ as a template codebase to be refactored and modularized, not a library to be imported.
- There is no unified CI/test harness in this repo yet.
  - typical test command pattern: nf-test test
- Use conda env files in `benchmark_project_script/conda_env/`:
  - `benchmark_py.yml` for Python method/evaluation scripts
  - `benchmark_R.yml` for R method/evaluation scripts
  - `eval_scib.yml` for `scib_metric_running.py`
- Typical command patterns:
  - Python integration scripts: `python <script>.py --all_targets <task...>`
  - R integration scripts: `Rscript <script>.R <task...>`
  - scIB evaluation: `python scib_metric_running.py --method <method> --target <task>`
- Prefer running commands from the directory layout expected by the scripts' relative paths.

## Conventions
- Task IDs are commonly hardcoded and include trailing underscores in integration scripts (for example `task4_`, `task6_`, `task8_`, `task9_`).
- Cross-species one-to-one ortholog mappings are task-specific and read from `../OrthoFinder/one2one_orthologs/*.csv` in many scripts.
- Common filtering defaults are repeated across methods: `min_genes=200`, `min_cells=10`.
- Metadata fields are normalized to `batch` and `celltype`; preserve existing label-handling logic when refactoring.
- Keep output locations consistent with existing conventions under `../output/method_outputs/` and `../output/evaluation/`.

## Pitfalls
- Many scripts depend on fragile relative paths and implicit working directories; avoid changing path semantics unless requested.
- R scripts include `setwd('./script')`, which assumes a specific folder layout that may differ from this repo snapshot. Validate cwd assumptions before modifying run instructions.
- Python and R flows use different file formats (`.h5ad` vs `.rds`); avoid introducing format-mismatch regressions.

## Documentation
- Original workflow and benchmark execution details: `benchmark_project_script/README.md`, use as a reference but not not copy explicitly, rather create new clean codebase and documentation in this repo.
- Repository-level direction and intended stack: `README.md`
- Follow link-first updates: extend existing docs rather than duplicating long procedural content in instruction files.
