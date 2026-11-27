nextflow.enable.types = true

/*
 * k-nearest-neighbour graph of an embedding into obsp/connectivities, obsp/distances and uns/neighbors: scanpy on
 * CPU, rapids-singlecell when the task has an accelerator. ext.args: --n-neighbors <int> (default 15), --n-dims
 * <int> (first columns only), --metric (default euclidean), --method umap|gauss, --seed (default 0). The GPU
 * search is exact; scanpy's is approximate on large data.
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SCANPY_NEIGHBORS {
    tag id
    label 'process_gpu'
    container task.accelerator ? 'quay.io/cellgeni/rapids-singlecell:0.18.1' : 'quay.io/cellgeni/scanpy:1.12.6'

    input:
    record(id: String, adata: Path, embedding: Path?, embedding_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, neighbors: file("${output_prefix ?: id}.neighbors.*"))

    topic:
    record(
        process: task.process,
        tool: task.accelerator ? 'rapids-singlecell' : 'scanpy',
        version: eval("scio.py version ${task.accelerator ? 'rapids-singlecell-cu12' : 'scanpy'}")
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    scanpy_neighbors.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${embedding ? "--embedding ${embedding}" : ''} \\
        ${embedding_key ? "--embedding-key ${embedding_key}" : ''} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.neighbors.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
