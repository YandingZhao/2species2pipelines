process FASTMNN_INTEGRATION {
    publishDir "${params.outdir}/fastmnn", mode: 'copy'

    input:
    tuple val(sample_id), path(input_a), path(input_b), val(species_a), val(species_b), path(fastmnn_script)

    output:
    path "${sample_id}_fastmnn_report.txt", emit: report
    path "${sample_id}_fastmnn_embedding.tsv", emit: embedding
    path "${sample_id}_fastmnn_integration.rds", emit: integrated_rds

    script:
    """
    Rscript ${fastmnn_script} \
      --input_a ${input_a} \
      --input_b ${input_b} \
      --sample_id ${sample_id} \
      --species_a ${species_a} \
      --species_b ${species_b}
    """

    stub:
    """
    printf "sample: ${sample_id}\n" > ${sample_id}_fastmnn_report.txt
    printf "species_a: ${species_a}\n" >> ${sample_id}_fastmnn_report.txt
    printf "species_b: ${species_b}\n" >> ${sample_id}_fastmnn_report.txt
    printf "status: stub_run\n" >> ${sample_id}_fastmnn_report.txt

    printf "cell\tfastmnn_1\tfastmnn_2\n" > ${sample_id}_fastmnn_embedding.tsv
    printf "stub_cell_1\t0.0\t0.0\n" >> ${sample_id}_fastmnn_embedding.tsv
    printf "stub fastmnn integration placeholder\n" > ${sample_id}_fastmnn_integration.rds
    """
}
