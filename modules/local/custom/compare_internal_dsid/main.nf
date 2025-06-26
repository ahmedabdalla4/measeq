/*
Exact matching of the N450 sequence to input DSId fasta file to determine if we already have the ID
    User has to provide the DSId fasta organized as:
    >DSId
    GTCAGTTCCA...

    Output is just a TSV that has:
    sample   matched_dsid   completeness
*/
process COMPARE_INTERNAL_DSID {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    path n450_fasta
    path id_fasta

    output:
    path "dsids.tsv", emit: dsid_tsv
    path "novel_dsids.tsv", optional: true, emit: novel_dsid_tsv
    path "versions.yml", emit: versions

    script:
    """
    compare_dsid.py \\
        --fasta $n450_fasta \\
        --dsid_fasta $id_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compare_dsid: x.x
    END_VERSIONS
    """

    stub:
    """
    touch dsids.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compare_dsid: x.x
    END_VERSIONS
    """
}
