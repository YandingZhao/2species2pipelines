process SEURAT4_INTEGRATION {
    publishDir "${params.outdir}/seurat4", mode: 'copy'

    input:
    tuple val(sample_id), path(input_a), path(input_b), val(species_a), val(species_b), path(seurat4_script)

    output:
    path "${sample_id}_seurat4_report.txt", emit: report
    path "${sample_id}_seurat4_embedding.tsv", emit: embedding
    path "${sample_id}_seurat4_integration.rds", emit: integrated_rds

    script:
    """
    Rscript ${seurat4_script} \
      --input_a ${input_a} \
      --input_b ${input_b} \
      --sample_id ${sample_id} \
      --species_a ${species_a} \
      --species_b ${species_b}
    """

    stub:
    """
    printf "sample: ${sample_id}\n" > ${sample_id}_seurat4_report.txt
    printf "species_a: ${species_a}\n" >> ${sample_id}_seurat4_report.txt
    printf "species_b: ${species_b}\n" >> ${sample_id}_seurat4_report.txt
    printf "status: stub_run\n" >> ${sample_id}_seurat4_report.txt

    printf "cell\tseurat4_1\tseurat4_2\n" > ${sample_id}_seurat4_embedding.tsv
    printf "stub_cell_1\t0.0\t0.0\n" >> ${sample_id}_seurat4_embedding.tsv
    printf "stub seurat4 integration placeholder\n" > ${sample_id}_seurat4_integration.rds
    """
}
