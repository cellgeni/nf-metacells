nextflow.enable.types = true

/*
 * Attach cell metadata to an AnnData store's obs and/or make its barcodes unique per sample, with adata-cli and
 * duckdb. Options go in ext.args: --barcode-column, --join inner|left, --delimiter, --sample-column,
 * --allow-dtype-change. Output format follows ext.extension: 'h5ad' (default) or 'zarr'.
 */
process ADATA_ATTACHOBS {
    tag id
    label 'process_single'
    container 'quay.io/cellgeni/adata-cli:0.17.6'

    input:
    record(id: String, adata: Path, metadata: Path?, sample_id: String?, output_prefix: String?, meta: Record?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(
        id: id,
        meta: meta,
        adata: file("${output_prefix ?: id}.${task.ext.extension ?: 'h5ad'}")
    )

    topic:
    record(process: task.process, tool: 'adata-cli', version: eval('adata --version')) >> 'versions'
    record(process: task.process, tool: 'duckdb', version: eval("duckdb --version | cut -d' ' -f1")) >> 'versions'

    script:
    def prefix = output_prefix ?: id
    def extension = task.ext.extension ?: 'h5ad'
    """
    adata_attachobs.sh \\
        --adata ${adata} \\
        --output ${prefix}.${extension} \\
        --sample-id ${sample_id ?: id} \\
        ${metadata ? "--metadata ${metadata}" : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def prefix = output_prefix ?: id
    def extension = task.ext.extension ?: 'h5ad'
    """
    ${extension == 'zarr' ? "mkdir ${prefix}.zarr" : "touch ${prefix}.${extension}"}
    """
}
