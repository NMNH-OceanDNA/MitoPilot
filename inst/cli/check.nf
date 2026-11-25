nextflow.enable.dsl=2

process check_mitopilot {

    executor params.curate.executor
    container params.curate.container
    clusterOptions {
        def opts = [
            (params.curate.executor == 'sge') ? '-S /bin/bash' : '',
            (params.curate.clusterOptions instanceof String) ? params.curate.clusterOptions : ''
        ].findAll { it }.join(' ')
        opts ?: null
    }

    tag "check"

    input:
    val(opts)

    output:
    stdout

    script:
    '''
    Rscript -e 'library(MitoPilot); cat("MitoPilot", as.character(packageVersion("MitoPilot")), "loads\\n")'
    '''
}

workflow {
    check_mitopilot(channel.of([cpus: 1, memory: 4]))
}
