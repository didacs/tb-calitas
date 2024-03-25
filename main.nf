nextflow.enable.dsl=2

// default parameters
params.spacers = "${workflow.projectDir}/test/spacers.txt"
params.fasta = "${workflow.projectDir}/test/chr21.sample.fa"
params.pam = 'nrg'
params.max_guide_diffs = 4
params.max_pam_mismatches = 0
params.max_gaps_between_guide_and_pam = 0
params.variants = null
params.gtf = "${workflow.projectDir}/test/chr21.sample.gtf"

//  processes
process get_spacers_metadata {
    publishDir "metadata", mode: 'copy'

    input:
        path spacers

    output:
        path "spacers_metadata.csv"

    script:
        """
        SP-benchling-query.py \
            --spacers ${spacers} \
            --outfile spacers_metadata.csv
        """
}


process run_calitas {
    input:
        tuple val(sp_id), val(bases)
    output:
        tuple val(sp_id), path("${sp_id}.txt")
    script:
        def variants_opt = ""
        if (params.variants != null) {
            variants_opt = "--variants ${params.variants}"
        }
        """
        calitas SearchReference \
            --guide ${bases}${params.pam} \
            --guide-id ${sp_id} \
            --ref ${params.fasta} \
            --max-guide-diffs ${params.max_guide_diffs} \
            --max-pam-mismatches ${params.max_pam_mismatches} \
            --max-gaps-between-guide-and-pam ${params.max_gaps_between_guide_and_pam} \
            --output ${sp_id}.txt \
            ${variants_opt}
        """
}

process annotate_sites {
    debug true
    publishDir "output", mode: 'copy'
    input:
        tuple val(sp_id), path(sites)
    output:
        path "${sp_id}.csv"
    script:
        """
        annotate_sites.py \
            --sites ${sites} \
            --gtf ${params.gtf} \
            --output ${sp_id}.csv
        """
}

workflow {
// get spacers metadata
    spacers_ch = Channel.fromPath(params.spacers)
    spacer_metadata = get_spacers_metadata(spacers_ch)
    spacer_metadata
        .splitCsv(header: true, sep: ',')
        .map(row -> tuple(row.sp_id, row.bases))
        .set {spacer_metadata_ch}
// run calitas
    sites_ch = run_calitas(spacer_metadata_ch)
// annotate sites
    annotate_sites(sites_ch)
}

workflow.onComplete {
    println "Pipeline execution summary:"
    println "Completed at: ${workflow.complete}"
    println "Duration    : ${workflow.duration}"
    println "Success     : ${workflow.success}"
    println "workDir     : ${workflow.workDir}"
    println "exit status : ${workflow.exitStatus}"
}