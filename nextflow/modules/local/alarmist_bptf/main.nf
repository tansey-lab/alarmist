process ALARMIST_BPTF {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/tansey-lab/alarmist:latest' :
        'ghcr.io/tansey-lab/alarmist:latest' }"

    input:
    tuple val(meta), path(patchify_results)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_bptf"
    """
    alarmist-bptf \\
        --input-dir ${patchify_results} \\
        --output-dir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: \$(python -c "import alarmist; print(alarmist.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_bptf"
    """
    mkdir -p ${prefix}
    touch ${prefix}/patch_loadings.parquet
    touch ${prefix}/lri_factors.parquet
    touch ${prefix}/model.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: 0.1.0
    END_VERSIONS
    """
}
