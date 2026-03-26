process ALARMIST_VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jeffquinnmsk/alarmist:latest' :
        'docker.io/jeffquinnmsk/alarmist:latest' }"

    input:
    tuple val(meta), path(glm_results), path(bptf_results), path(project_results), path(patchify_results)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_visualize"
    """
    alarmist-visualize \\
        --glm-dir ${glm_results} \\
        --bptf-dir ${bptf_results} \\
        --project-dir ${project_results} \\
        --patchify-dir ${patchify_results} \\
        --output-dir ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: \$(python -c "import alarmist; print(alarmist.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_visualize"
    """
    mkdir -p ${prefix}
    touch ${prefix}/volcano_motif_0_mqc.png
    touch ${prefix}/forest_motif_0_mqc.png
    touch ${prefix}/motif_heatmap_mqc.png
    touch ${prefix}/motif_distributions_mqc.png
    touch ${prefix}/bptf_diagnostics_mqc.png
    touch ${prefix}/spatial_celltypes_mqc.png
    touch ${prefix}/spatial_motif_0_loading_mqc.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: 0.1.0
    END_VERSIONS
    """
}
