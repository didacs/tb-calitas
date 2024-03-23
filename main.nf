nextflow.enable.dsl=2

// default parameters
params.spacers = "${workflow.projectDir}/test/spacers.txt"
params.pam = 'nrg'
params.fasta = "${workflow.projectDir}/test/chr21.sample.fa"
params.assembly_name = "chr21"
params.species = "human"
params.variants = null
params.max_guide_diffs = 4
params.max_pam_mismatches = 0
params.max_gaps_between_guide_and_pam = 0

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


process prepare_ref {
    input:
        path fasta
        val assembly_name
        val species
    output:
        path fasta, emit: 'fasta'
        path "${fasta}.fai", emit: 'index'
        path "${fasta}.dict", emit: 'dict'
    script:
    """
    samtools faidx ${fasta} && \
    samtools dict -a ${assembly_name} -s ${species} -o ${fasta}.dict ${fasta}
    """
}


process run_calitas {
    publishDir "output", mode: 'copy'
    input:
        tuple val(sp_id), val(bases)
        path fasta
        path index
        path dict
    output:
        path "${sp_id}.calitas.txt"
    script:
        def variants_opt = ""
        if (params.variants != null) {
            variants_opt = "--variants ${params.variants}"
        }
        """
        calitas SearchReference \
            --guide ${bases}${params.pam} \
            --guide-id ${sp_id} \
            --ref ${fasta} \
            --max-guide-diffs ${params.max_guide_diffs} \
            --max-pam-mismatches ${params.max_pam_mismatches} \
            --max-gaps-between-guide-and-pam ${params.max_gaps_between_guide_and_pam} \
            --output ${sp_id}.calitas.txt \
            ${variants_opt}
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
// prepare reference
    fasta_ch = Channel.fromPath(params.fasta)
    fasta_ch.view()
    ref_ch = prepare_ref(fasta_ch, params.assembly_name, params.species)
// run calitas
    run_calitas(spacer_metadata_ch, ref_ch.fasta.first(), ref_ch.index.first(), ref_ch.dict.first())
}

workflow.onComplete {
    println "Pipeline execution summary:"
    println "Completed at: ${workflow.complete}"
    println "Duration    : ${workflow.duration}"
    println "Success     : ${workflow.success}"
    println "workDir     : ${workflow.workDir}"
    println "exit status : ${workflow.exitStatus}"
}