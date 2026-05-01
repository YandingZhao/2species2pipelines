nextflow.enable.dsl = 2

include { NFCORE_BASE } from './workflows/nfcore_base'

workflow {
    NFCORE_BASE()
}
