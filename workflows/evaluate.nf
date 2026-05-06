nextflow.enable.dsl = 2

include { EVALUATE_INTEGRATION } from '../modules/local/evaluate_integration'
include { AGGREGATE_UNSCALED_METRICS } from '../modules/local/aggregate_unscaled_metrics'

workflow EVALUATE {
    take:
    ch_integrated_h5ad

    main:
    evaluate_integration_script = file("${projectDir}/scripts/run_evaluate_integration.py")
    aggregate_unscaled_metrics_script = file("${projectDir}/scripts/aggregate_unscaled_metrics.py")

    ch_h5ad_eval = ch_integrated_h5ad.map { integrated_h5ad ->
        tuple(integrated_h5ad, evaluate_integration_script)
    }

    EVALUATE_INTEGRATION(ch_h5ad_eval)

    ch_metrics_to_aggregate = EVALUATE_INTEGRATION.out.metrics.collect()
    ch_aggregation_script = channel.value(aggregate_unscaled_metrics_script)
    AGGREGATE_UNSCALED_METRICS(ch_metrics_to_aggregate, ch_aggregation_script)

    emit:
    evaluation_reports = EVALUATE_INTEGRATION.out.report
    evaluation_metrics = EVALUATE_INTEGRATION.out.metrics
    evaluation_metrics_scaled = EVALUATE_INTEGRATION.out.metrics_scaled
    evaluation_metrics_combined_report = AGGREGATE_UNSCALED_METRICS.out.report
    evaluation_metrics_combined_long = AGGREGATE_UNSCALED_METRICS.out.long
    evaluation_metrics_combined_image = AGGREGATE_UNSCALED_METRICS.out.image
}
