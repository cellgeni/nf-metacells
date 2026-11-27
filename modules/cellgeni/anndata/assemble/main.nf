nextflow.enable.types = true

/*
 * Assemble a base AnnData object and step outputs (raw directories, slim or full objects) into one object: a copy
 * of the base with every element of the sources written into it (later sources win; elements the base already has
 * are kept from the base for object sources). The step manifests go to uns['nf_metacells'].
 * Output format: ext.adata_format, h5ad (default) or zarr.
 */
process ANNDATA_ASSEMBLE {
    tag id
    label 'process_single'
    container 'quay.io/cellgeni/scanpy:1.12.6'

    input:
    record(id: String, adata: Path, sources: List<Path>, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, assembled: file("${output_prefix ?: id}.assemble.*"))

    topic:
    record(process: task.process, tool: 'anndata', version: eval('scio.py version anndata')) >> 'versions'

    script:
    """
    anndata_assemble.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${sources ? "--sources ${sources.join(' ')}" : ''} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'}
    """

    stub:
    def format = task.ext.adata_format ?: 'h5ad'
    """
    ${format == 'zarr' ? "mkdir ${output_prefix ?: id}.assemble.zarr" : "touch ${output_prefix ?: id}.assemble.h5ad"}
    """
}
