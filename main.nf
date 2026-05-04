nextflow.enable.dsl = 2

include { INTEGRATE } from './workflows/integrate'

workflow {
    INTEGRATE()
}
