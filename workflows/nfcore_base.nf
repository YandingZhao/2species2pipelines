nextflow.enable.dsl = 2

include { MAKE_RUN_METADATA } from '../modules/local/make_run_metadata'
include { HARMONY_INTEGRATION } from '../modules/local/harmony_integration'

workflow NFCORE_BASE {
    main:
    ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    ch_checked = ch_samples.ifEmpty { error "No rows found in input samplesheet: ${params.input}" }
    ch_reports = ch_checked.map { row -> tuple(row.sample as String, row) }
    harmony_script = file("${projectDir}/scripts/run_harmony_module.R")
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

    MAKE_RUN_METADATA(ch_reports)
    HARMONY_INTEGRATION(ch_harmony)

    emit:
    report_files = MAKE_RUN_METADATA.out.report
    harmony_reports = HARMONY_INTEGRATION.out.report
    harmony_pca = HARMONY_INTEGRATION.out.pca
    harmony_rds = HARMONY_INTEGRATION.out.integrated_rds
}
