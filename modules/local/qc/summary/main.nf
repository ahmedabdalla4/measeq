// Summary QC
//  Using artic env for the moment as it has a bunch of tools and
//  is used in earlier steps
process MAKE_FINAL_QC_CSV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9941f0d0ff90a15b0f7804a7487c30957fb6e6cf:aa1fff77714712ab105e19c770efcb635ed300ff-0' :
        'biocontainers/mulled-v2-9941f0d0ff90a15b0f7804a7487c30957fb6e6cf:aa1fff77714712ab105e19c770efcb635ed300ff-0' }"

    input:
    path concat_csv
    path metadata
    val neg_control_pct_threshold
    val neg_ctrl_substrings
    val skip_negative_grading

    output:
    path "overall.qc.csv", emit: csv
    path "overall.xlsx", emit: xlsx
    path "versions.yml", emit: versions

    script:
    def version = workflow.manifest.version
    def metadataArg = metadata ? "--metadata $metadata" : ""
    def negativeGradingArg = skip_negative_grading ? "--skip_ctrl_grading" : ""
    """
    summary_qc.py \\
        --csv $concat_csv \\
        $metadataArg \\
        $negativeGradingArg \\
        --threshold $neg_control_pct_threshold \\
        --neg_ctrl_substrings '$neg_ctrl_substrings' \\
        --version $version

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch overall.qc.csv
    touch overall.xlsx

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
