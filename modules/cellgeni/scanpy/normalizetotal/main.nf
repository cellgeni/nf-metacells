nextflow.enable.types = true

/*
 * Normalise counts to the same total per cell into layers/normalized: scanpy on CPU, rapids-singlecell when the
 * task has an accelerator. Counts come from counts_key of the base object (X, layers/<name> or raw/X; default X).
 * ext.args: --target-sum, --exclude-highly-expressed, --max-fraction.
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SCANPY_NORMALIZETOTAL {
    tag id
    label 'process_gpu'
    container task.accelerator ?
        'quay.io/cellgeni/rapids-singlecell:0.18.1' :
        'community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90'

    input:
    record(id: String, adata: Path, counts_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, normalized: file("${output_prefix ?: id}.normalizetotal.*"))

    topic:
    record(
        process: task.process,
        tool: task.accelerator ? 'rapids-singlecell' : 'scanpy',
        version: eval("scio.py version ${task.accelerator ? 'rapids-singlecell-cu12' : 'scanpy'}")
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    export PYTHONNOUSERSITE=1 NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}/numba MPLCONFIGDIR=\${TMPDIR:-/tmp}/matplotlib
    scanpy_normalizetotal.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        --counts-key ${counts_key ?: 'X'} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.normalizetotal.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
