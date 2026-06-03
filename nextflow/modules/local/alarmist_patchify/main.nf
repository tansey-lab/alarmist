process ALARMIST_PATCHIFY {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jeffquinnmsk/alarmist:latest' :
        'docker.io/jeffquinnmsk/alarmist:latest' }"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_args = params.sample_column ? "--multi-sample --sample-column ${params.sample_column}" : ''
    prefix = task.ext.prefix ?: "${meta.id}_patchify"
    """
    alarmist-patchify \\
        --adata ${adata} \\
        --output-dir ${prefix} \\
        --cell-type-column ${params.cell_type_column} \\
        --patch-size ${params.patch_size} \\
        --resource ${params.resource} \\
        ${sample_args} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: \$(python -c "import alarmist; print(alarmist.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_patchify"
    """
    mkdir -p ${prefix}
    touch ${prefix}/patch_lri_matrix.parquet
    touch ${prefix}/patch_metadata.parquet
    touch ${prefix}/lri_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alarmist: 0.1.0
    END_VERSIONS
    """
}
