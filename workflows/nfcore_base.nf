nextflow.enable.dsl = 2

include { MAKE_RUN_METADATA } from '../modules/local/make_run_metadata'
include { HARMONY_INTEGRATION } from '../modules/local/harmony_integration'
include { SEURAT4_INTEGRATION } from '../modules/local/seurat4_integration'
include { FASTMNN_INTEGRATION } from '../modules/local/fastmnn_integration'
include { BBKNN_INTEGRATION } from '../modules/local/bbknn_integration'

workflow NFCORE_BASE {
    main:
    ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    ch_checked = ch_samples.ifEmpty { error "No rows found in input samplesheet: ${params.input}" }
    ch_reports = ch_checked.map { row -> tuple(row.sample as String, row) }
    ch_bbknn_samples = channel
        .fromPath(params.bbknn_input, checkIfExists: true)
        .splitCsv(header: true)
        .ifEmpty { error "No rows found in BBKNN samplesheet: ${params.bbknn_input}" }

    harmony_script = file("${projectDir}/scripts/run_harmony_module.R")
    seurat4_script = file("${projectDir}/scripts/run_seurat4_module.R")
    fastmnn_script = file("${projectDir}/scripts/run_fastmnn_module.R")
    bbknn_script = file("${projectDir}/scripts/run_bbknn_module.py")
    ch_harmony = ch_checked.map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            harmony_script
        )
    }
    ch_seurat4 = ch_checked.map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            seurat4_script
        )
    }
    ch_fastmnn = ch_checked.map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            fastmnn_script
        )
    }
    ch_bbknn = ch_bbknn_samples.map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            bbknn_script
        )
    }

    MAKE_RUN_METADATA(ch_reports)
    HARMONY_INTEGRATION(ch_harmony)
    SEURAT4_INTEGRATION(ch_seurat4)
    FASTMNN_INTEGRATION(ch_fastmnn)
    BBKNN_INTEGRATION(ch_bbknn)

    emit:
    report_files = MAKE_RUN_METADATA.out.report
    harmony_reports = HARMONY_INTEGRATION.out.report
    harmony_pca = HARMONY_INTEGRATION.out.pca
    harmony_rds = HARMONY_INTEGRATION.out.integrated_rds
    seurat4_reports = SEURAT4_INTEGRATION.out.report
    seurat4_pca = SEURAT4_INTEGRATION.out.pca
    seurat4_rds = SEURAT4_INTEGRATION.out.integrated_rds
    fastmnn_reports = FASTMNN_INTEGRATION.out.report
    fastmnn_embedding = FASTMNN_INTEGRATION.out.embedding
    fastmnn_rds = FASTMNN_INTEGRATION.out.integrated_rds
    bbknn_reports = BBKNN_INTEGRATION.out.report
    bbknn_pca = BBKNN_INTEGRATION.out.pca
    bbknn_h5ad = BBKNN_INTEGRATION.out.integrated_h5ad
}
