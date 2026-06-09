process ALARMIST_GLM {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jeffquinnmsk/alarmist:latest' :
        'docker.io/jeffquinnmsk/alarmist:latest' }"

    input:
    tuple val(meta), path(project_results), path(patchify_results)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_glm"
    """
    alarmist-glm \\
        --input-dir ${project_results} \\
        --adata ${project_results}/projected_adata.h5ad \\
        --patch-lri-dir ${patchify_results} \\
        --output-dir ${prefix} \\
        --cell-type-column ${params.cell_type_column} \\
        ${params.prefilter_spearman ? '--prefilter-spearman' : '--no-prefilter-spearman'} \\
        --spearman-pval-threshold ${params.spearman_pval_threshold} \\
        --spearman-chunk-size ${params.spearman_chunk_size} \\
        --backend ${params.glm_backend} \\
        --device ${params.glm_device} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: \$(python -c "import alarmist; print(alarmist.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_glm"
    """
    mkdir -p ${prefix}
    touch ${prefix}/glm_results.parquet
    touch ${prefix}/coefficients.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: 0.1.0
    END_VERSIONS
    """
}
