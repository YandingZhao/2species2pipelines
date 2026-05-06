process AGGREGATE_UNSCALED_METRICS {
        publishDir "${params.outdir}/evaluation", mode: 'copy'
        tag "aggregate_unscaled_metrics"

        input:
        path metrics_files
        path aggregation_script

        output:
        path "combined_unscaled_metrics_report.tsv", emit: report
        path "combined_unscaled_metrics_long.tsv", emit: long
        path "combined_unscaled_metrics_report.png", emit: image

        script:
        def metric_args = metrics_files.collect { metric_file -> "${metric_file}" }.join(' ')
        """
        python ${aggregation_script} \
            --metrics_files ${metric_args} \
            --output_report combined_unscaled_metrics_report.tsv \
            --output_long combined_unscaled_metrics_long.tsv \
            --output_figure combined_unscaled_metrics_report.png
        """

        stub:
        """
        printf "integration\tembedding\tstatus\treason\n" > combined_unscaled_metrics_report.tsv
        printf "stub\tNA\tstub_run\tstub_run\n" >> combined_unscaled_metrics_report.tsv
        printf "integration\tembedding\tmetric\tvalue\n" > combined_unscaled_metrics_long.tsv
        printf "stub\tNA\tstatus\tstub_run\n" >> combined_unscaled_metrics_long.tsv
        touch combined_unscaled_metrics_report.png
        """
}
