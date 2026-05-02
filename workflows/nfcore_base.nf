nextflow.enable.dsl = 2

include { MAKE_RUN_METADATA } from '../modules/local/make_run_metadata'
include { HARMONY_INTEGRATION } from '../modules/local/harmony_integration'
include { SEURAT4_INTEGRATION } from '../modules/local/seurat4_integration'
include { FASTMNN_INTEGRATION } from '../modules/local/fastmnn_integration'
include { BBKNN_INTEGRATION } from '../modules/local/bbknn_integration'
include { SCANORAMA_INTEGRATION } from '../modules/local/scanorama_integration'
include { SCVI_INTEGRATION } from '../modules/local/scvi_integration'
include { SEURAT_TO_ANNDATA_PAIR } from '../modules/local/seurat_to_anndata_pair'
include { ORTHOLOG_CONVERT_PAIR } from '../modules/local/ortholog_convert_pair'

workflow NFCORE_BASE {
    main:
    ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    ch_checked = ch_samples.ifEmpty { error "No rows found in input samplesheet: ${params.input}" }
    ch_reports = ch_checked.map { row -> tuple(row.sample as String, row) }

    harmony_script = file("${projectDir}/scripts/run_harmony_module.R")
    seurat4_script = file("${projectDir}/scripts/run_seurat4_module.R")
    fastmnn_script = file("${projectDir}/scripts/run_fastmnn_module.R")
    bbknn_script = file("${projectDir}/scripts/run_bbknn_module.py")
    scanorama_script = file("${projectDir}/scripts/run_scanorama_module.py")
    scvi_script = file("${projectDir}/scripts/run_scvi_module.py")
    seurat_to_anndata_script = file("${projectDir}/scripts/run_seurat_to_anndata_pair.R")
    ortholog_convert_script = file("${projectDir}/scripts/run_ortholog_convert_pair.R")

    ch_rds_for_ortholog = ch_checked
        .filter { row ->
            row.source_a.toString().toLowerCase().endsWith('.rds') &&
            row.source_b.toString().toLowerCase().endsWith('.rds')
        }
        .map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            ortholog_convert_script
        )
    }

    ORTHOLOG_CONVERT_PAIR(ch_rds_for_ortholog)

    ch_rds_converted = ORTHOLOG_CONVERT_PAIR.out.rds_pair

    ch_harmony = ch_rds_converted.map { sample_id, input_a, input_b, species_a, species_b ->
        tuple(
            sample_id,
            input_a,
            input_b,
            species_a,
            species_b,
            harmony_script
        )
    }
    ch_seurat4 = ch_rds_converted.map { sample_id, input_a, input_b, species_a, species_b ->
        tuple(
            sample_id,
            input_a,
            input_b,
            species_a,
            species_b,
            seurat4_script
        )
    }
    ch_fastmnn = ch_rds_converted.map { sample_id, input_a, input_b, species_a, species_b ->
        tuple(
            sample_id,
            input_a,
            input_b,
            species_a,
            species_b,
            fastmnn_script
        )
    }
    ch_bbknn_h5ad = ch_checked
        .filter { row ->
            row.source_a.toString().toLowerCase().endsWith('.h5ad') &&
            row.source_b.toString().toLowerCase().endsWith('.h5ad')
        }
        .map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            bbknn_script
        )
    }

    ch_bbknn_rds_convert = ch_rds_converted.map { sample_id, input_a, input_b, species_a, species_b ->
        tuple(
            sample_id,
            input_a,
            input_b,
            species_a,
            species_b,
            seurat_to_anndata_script
        )
    }

    SEURAT_TO_ANNDATA_PAIR(ch_bbknn_rds_convert)

    ch_bbknn_from_rds = SEURAT_TO_ANNDATA_PAIR.out.bbknn_input.map { sample_id, input_a_h5ad, input_b_h5ad, species_a, species_b ->
        tuple(sample_id, input_a_h5ad, input_b_h5ad, species_a, species_b, bbknn_script)
    }

    ch_bbknn = ch_bbknn_h5ad.mix(ch_bbknn_from_rds)
        .ifEmpty { error "No valid BBKNN inputs found. Provide both inputs as .h5ad or both as .rds in ${params.input}" }

    ch_scanorama_h5ad = ch_checked
        .filter { row ->
            row.source_a.toString().toLowerCase().endsWith('.h5ad') &&
            row.source_b.toString().toLowerCase().endsWith('.h5ad')
        }
        .map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            scanorama_script
        )
    }

    ch_scanorama_from_rds = SEURAT_TO_ANNDATA_PAIR.out.bbknn_input.map { sample_id, input_a_h5ad, input_b_h5ad, species_a, species_b ->
        tuple(sample_id, input_a_h5ad, input_b_h5ad, species_a, species_b, scanorama_script)
    }

    ch_scanorama = ch_scanorama_h5ad.mix(ch_scanorama_from_rds)
        .ifEmpty { error "No valid Scanorama inputs found. Provide both inputs as .h5ad or both as .rds in ${params.input}" }

    ch_scvi_h5ad = ch_checked
        .filter { row ->
            row.source_a.toString().toLowerCase().endsWith('.h5ad') &&
            row.source_b.toString().toLowerCase().endsWith('.h5ad')
        }
        .map { row ->
        tuple(
            row.sample as String,
            file(row.source_a),
            file(row.source_b),
            (row.species_a ?: "unknown") as String,
            (row.species_b ?: "unknown") as String,
            scvi_script
        )
    }

    ch_scvi_from_rds = SEURAT_TO_ANNDATA_PAIR.out.bbknn_input.map { sample_id, input_a_h5ad, input_b_h5ad, species_a, species_b ->
        tuple(sample_id, input_a_h5ad, input_b_h5ad, species_a, species_b, scvi_script)
    }

    ch_scvi = ch_scvi_h5ad.mix(ch_scvi_from_rds)
        .ifEmpty { error "No valid scVI inputs found. Provide both inputs as .h5ad or both as .rds in ${params.input}" }

    MAKE_RUN_METADATA(ch_reports)
    HARMONY_INTEGRATION(ch_harmony)
    SEURAT4_INTEGRATION(ch_seurat4)
    FASTMNN_INTEGRATION(ch_fastmnn)
    BBKNN_INTEGRATION(ch_bbknn)
    SCANORAMA_INTEGRATION(ch_scanorama)
    SCVI_INTEGRATION(ch_scvi)

    emit:
    report_files = MAKE_RUN_METADATA.out.report
    ortholog_reports = ORTHOLOG_CONVERT_PAIR.out.report
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
    scanorama_reports = SCANORAMA_INTEGRATION.out.report
    scanorama_pca = SCANORAMA_INTEGRATION.out.pca
    scanorama_h5ad = SCANORAMA_INTEGRATION.out.integrated_h5ad
    scvi_reports = SCVI_INTEGRATION.out.report
    scvi_pca = SCVI_INTEGRATION.out.pca
    scvi_h5ad = SCVI_INTEGRATION.out.integrated_h5ad
}
