process ORTHOLOG_CONVERT_PAIR {
    publishDir "${params.outdir}/ortholog_conversion", mode: 'copy'

    input:
    tuple val(sample_id), path(input_a), path(input_b), val(species_a), val(species_b), path(converter_script)

    output:
    tuple val(sample_id), path("${sample_id}_a_ortholog.rds"), path("${sample_id}_b_ortholog.rds"), val(species_a), val(species_b), emit: rds_pair
    path "${sample_id}_ortholog_report.txt", emit: report

    script:
    """
    Rscript ${converter_script} \
      --input_a ${input_a} \
      --input_b ${input_b} \
      --sample_id ${sample_id} \
      --species_a ${species_a} \
      --species_b ${species_b}
    """

    stub:
    """
    printf "stub ortholog a\n" > ${sample_id}_a_ortholog.rds
    printf "stub ortholog b\n" > ${sample_id}_b_ortholog.rds
    printf "sample: ${sample_id}\n" > ${sample_id}_ortholog_report.txt
    printf "status: stub_run\n" >> ${sample_id}_ortholog_report.txt
    """
}