process SCANORAMA_INTEGRATION {
    publishDir "${params.outdir}/scanorama", mode: 'copy'

    input:
    tuple val(sample_id), path(input_a), path(input_b), val(species_a), val(species_b), path(scanorama_script)

    output:
    path "${sample_id}_scanorama_report.txt", emit: report
    path "${sample_id}_scanorama_embedding.tsv", emit: embedding
    path "${sample_id}_scanorama_integration.h5ad", emit: integrated_h5ad

    script:
    """
    python ${scanorama_script} \
      --input_a ${input_a} \
      --input_b ${input_b} \
      --sample_id ${sample_id} \
      --species_a ${species_a} \
      --species_b ${species_b}
    """

    stub:
    """
    printf "sample: ${sample_id}\n" > ${sample_id}_scanorama_report.txt
    printf "species_a: ${species_a}\n" >> ${sample_id}_scanorama_report.txt
    printf "species_b: ${species_b}\n" >> ${sample_id}_scanorama_report.txt
    printf "status: stub_run\n" >> ${sample_id}_scanorama_report.txt

    printf "cell\tPC1\tPC2\n" > ${sample_id}_scanorama_embedding.tsv
    printf "stub_cell_1\t0.0\t0.0\n" >> ${sample_id}_scanorama_embedding.tsv

    printf "stub scanorama integration placeholder\n" > ${sample_id}_scanorama_integration.h5ad
    """
}