nextflow.enable.types = true

/*
 * Convert a single-cell matrix (10x .h5, 10x mtx folder, h5ad, zarr, loom, ...) to an AnnData store with
 * adata-cli. Output format follows ext.extension: 'h5ad' (default) or 'zarr'.
 */
process ADATA_CONVERT {
    tag id
    label 'process_single'
    container 'quay.io/cellgeni/adata-cli:0.17.6'

    input:
    record(id: String, source: Path, convert_args: String?, output_prefix: String?, meta: Record?)

    stage:
    stageAs source, 'input/*'

    output:
    record(
        id: id,
        meta: meta,
        adata: file("${output_prefix ?: id}.${task.ext.extension ?: 'h5ad'}")
    )

    topic:
    record(process: task.process, tool: 'adata-cli', version: eval('adata --version')) >> 'versions'

    script:
    def prefix = output_prefix ?: id
    def extension = task.ext.extension ?: 'h5ad'
    """
    adata convert ${source} -o ${prefix}.${extension} ${convert_args ?: ''} ${task.ext.args ?: ''}
    """

    stub:
    def prefix = output_prefix ?: id
    def extension = task.ext.extension ?: 'h5ad'
    """
    ${extension == 'zarr' ? "mkdir ${prefix}.zarr" : "touch ${prefix}.${extension}"}
    """
}
