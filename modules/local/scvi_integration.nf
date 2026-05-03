process SCVI_INTEGRATION {
    publishDir "${params.outdir}/scvi", mode: 'copy'

    input:
    tuple val(sample_id), path(input_a), path(input_b), val(species_a), val(species_b), path(scvi_script)

    output:
    path "${sample_id}_scvi_report.txt", emit: report
    path "${sample_id}_scvi_pca.tsv", emit: pca
    path "${sample_id}_scvi_integration.h5ad", emit: integrated_h5ad

    script:
    """
    python ${scvi_script} \
      --input_a ${input_a} \
      --input_b ${input_b} \
      --sample_id ${sample_id} \
      --species_a ${species_a} \
      --species_b ${species_b}
    """

    stub:
    """
    printf "sample: ${sample_id}\n" > ${sample_id}_scvi_report.txt
    printf "species_a: ${species_a}\n" >> ${sample_id}_scvi_report.txt
    printf "species_b: ${species_b}\n" >> ${sample_id}_scvi_report.txt
    printf "status: stub_run\n" >> ${sample_id}_scvi_report.txt

    printf "cell\tPC1\tPC2\n" > ${sample_id}_scvi_pca.tsv
    printf "stub_cell_1\t0.0\t0.0\n" >> ${sample_id}_scvi_pca.tsv

    printf "stub scvi integration placeholder\n" > ${sample_id}_scvi_integration.h5ad
    """
}