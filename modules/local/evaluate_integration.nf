process EVALUATE_INTEGRATION {
        publishDir "${params.outdir}/evaluation", mode: 'copy'
        tag "${integrated_h5ad.baseName}"

        input:
        tuple path(integrated_h5ad), path(evaluation_script)

        output:
        path "${integrated_h5ad.baseName}_scib_report.txt", emit: report
        path "${integrated_h5ad.baseName}_scib_metrics.tsv", emit: metrics

        script:
        """
        python ${evaluation_script} \
            --input_h5ad ${integrated_h5ad} \
            --batch_key batch \
            --label_key celltype
        """

        stub:
        """
        printf "status: stub_run\n" > ${integrated_h5ad.baseName}_scib_report.txt
        printf "metric\tvalue\n" > ${integrated_h5ad.baseName}_scib_metrics.tsv
        printf "status\tstub_run\n" >> ${integrated_h5ad.baseName}_scib_metrics.tsv
        """
}