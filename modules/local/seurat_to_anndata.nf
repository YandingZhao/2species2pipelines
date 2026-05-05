process SEURAT_TO_ANNDATA_PAIR {
    publishDir "${params.outdir}/seurat_conversion", mode: 'copy'

    input:
    tuple val(sample_id), path(input_a), path(input_b), val(species_a), val(species_b), path(converter_script)

    output:
    tuple val(sample_id), path("${input_a.baseName}.h5ad"), path("${input_b.baseName}.h5ad"), val(species_a), val(species_b), emit: anndata_pairs

    script:
    """
    Rscript ${converter_script} \
      --input_a ${input_a} \
      --input_b ${input_b} \
      --sample_id ${sample_id}
    """

    stub:
    """
    printf "stub h5ad from rds a\n" > ${input_a.baseName}.h5ad
    printf "stub h5ad from rds b\n" > ${input_b.baseName}.h5ad
    """
}

process SEURAT_TO_ANNDATA_SINGLE {
    publishDir "${params.outdir}/seurat_conversion", mode: 'copy'

    input:
    tuple path(input_rds), path(converter_script)

    output:
    path("${input_rds.baseName}.h5ad"), emit: anndata_single

    script:
    """
    Rscript ${converter_script} \
      --input_a ${input_rds} \
    """

    stub:
    """
    printf "stub h5ad from rds single\n" > ${input_rds.baseName}.h5ad
    """
}