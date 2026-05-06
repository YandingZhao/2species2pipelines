nextflow.enable.dsl = 2

include { INTEGRATE } from './workflows/integrate'
include { EVALUATE } from './workflows/evaluate'

workflow {
    INTEGRATE()
    EVALUATE(INTEGRATE.out.integrated_h5ad)
}
