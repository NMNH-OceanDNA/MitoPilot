nextflow.enable.dsl=2

// Curate one lifted annotation (plus its inert pad row) with the released
// MitoPilot core function; the app validates it. Launched by the app from a
// scratch folder under <project>/.liftover/; never part of a pipeline run.

params.lift_annotations = null
params.lift_assembly = null
params.lift_coverage = null
params.lift_ref_db = null
params.lift_blast_ref = null
params.lift_params = null
params.lift_genetic_code = 2
params.lift_max_blast_hits = 10
params.lift_outdir = null

process liftover_curate {

    executor 'local'
    errorStrategy 'terminate'
    publishDir "${params.lift_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(opts), path(annotations), path(assembly), path(coverage), path(ref_db), path(blast_ref_files, stageAs: 'blast_ref_*.json')

    output:
    path "curate"

    shell:
    cov = coverage.name == 'NO_FILE' ? 'NULL' : "'${coverage}'"
    refs = (blast_ref_files instanceof List ? blast_ref_files : [blast_ref_files]).join(' ')
    '''
    export OMP_NUM_THREADS=1
    mkdir -p curate
    # container engines may write as root; keep the scratch folder removable
    trap 'chmod -R a+rwX curate 2>/dev/null || true' EXIT
    Rscript -e "withCallingHandlers( \
        MitoPilot:::curate_mito_core( \
            annotations_fn = '!{annotations}', \
            assembly_fn = '!{assembly}', \
            coverage_fn = !{cov}, \
            genetic_code = !{params.lift_genetic_code}, \
            out_dir = 'curate', \
            max_blast_hits = !{params.lift_max_blast_hits}, \
            params = '!{params.lift_params}', \
            ref_dir = '!{ref_db}', \
            blast_ref_file = '!{refs}', \
            feature_trim = FALSE, \
            flip_rRNA_minus_strand = FALSE, \
            ref_based_rc = FALSE \
        ), \
        warning = function(w) if (grepl('had status', conditionMessage(w))) stop(conditionMessage(w), call. = FALSE) \
    )"
    '''
}

workflow {
    liftover_curate(channel.of(tuple(
        [cpus: 1, memory: 2],  // for project-config closures that read opts
        file(params.lift_annotations),
        file(params.lift_assembly),
        file(params.lift_coverage),
        file(params.lift_ref_db, type: 'dir'),
        params.lift_blast_ref.tokenize(',').collect { file(it) }
    )))
}
