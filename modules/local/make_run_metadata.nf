process MAKE_RUN_METADATA {
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample_id), val(row)

    output:
    path "${sample_id}.txt", emit: report

    script:
    """
    cat > ${sample_id}.txt <<'EOF'
    sample: ${sample_id}
    species_a: ${row.species_a}
    species_b: ${row.species_b}
    source_a: ${row.source_a}
    source_b: ${row.source_b}
    EOF
    """
}
