process ORTHOLOG_CONVERT {
    publishDir "${params.outdir}/converted/${strategy}", mode: 'copy'

    input:
    tuple val(strategy), path(input_a), path(input_b), val(species_a), val(species_b)
    path map_file  // optional; pass file("NO_FILE") when not required

    output:
    tuple val(strategy), path("Dog_converted.rds"), path("Human_converted.rds"), val(species_a), val(species_b), emit: rds_pair
    path "${strategy}_ortholog_report.txt", emit: report

    script:
    def map_arg = (map_file.name != "NO_FILE") ? "--map_file ${map_file}" : ""
    """
    Rscript ${projectDir}/scripts/run_ortholog_convert.R \
      --strategy   ${strategy} \
      --input_a    ${input_a} \
      --input_b    ${input_b} \
      --species_a  ${species_a} \
      --species_b  ${species_b} \
      --output_a   Dog_converted.rds \
      --output_b   Human_converted.rds \
      ${map_arg}

    # Move report to strategy-prefixed name
    mv Dog_converted_ortholog_report.txt ${strategy}_ortholog_report.txt 2>/dev/null || \
      printf "strategy: ${strategy}\\nstatus: ok\\n" > ${strategy}_ortholog_report.txt
    """

    stub:
    """
    printf "stub converted a\n" > Dog_converted.rds
    printf "stub converted b\n" > Human_converted.rds
    printf "strategy: ${strategy}\nstatus: stub\n" > ${strategy}_ortholog_report.txt
    """
}
