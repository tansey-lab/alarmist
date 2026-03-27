process ALARMIST_PROJECT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jeffquinnmsk/alarmist:latest' :
        'docker.io/jeffquinnmsk/alarmist:latest' }"

    input:
    tuple val(meta), path(adata), path(patchify_results), path(bptf_results)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_project"
    """
    alarmist-project \\
        --adata ${adata} \\
        --bptf-dir ${bptf_results} \\
        --patch-lri-dir ${patchify_results} \\
        --output-dir ${prefix} \\
        --cell-type-column ${params.cell_type_column} \\
        --resource ${params.resource} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: \$(python -c "import alarmist; print(alarmist.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_project"
    """
    mkdir -p ${prefix}
    touch ${prefix}/cell_motif_scores.parquet
    touch ${prefix}/projected_adata.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: 0.1.0
    END_VERSIONS
    """
}
