nextflow.enable.types = true

/*
 * Log-transforms a normalised layer, log(1 + x), into layers/log1p: scanpy on CPU, rapids-singlecell when the
 * task has an accelerator. ext.args: --layer-key <element> (default layers/normalized), --base <float> (default:
 * natural logarithm).
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SCANPY_LOG1P {
    tag id
    label 'process_gpu'
    container task.accelerator ?
        'quay.io/cellgeni/rapids-singlecell:0.18.1' :
        'community.wave.seqera.io/library/python_pip_scanpy_anndata_pruned:a5329e1ed934ad90'

    input:
    record(id: String, adata: Path, normalized: Path?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, log1p: file("${output_prefix ?: id}.log1p.*"))

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
    scanpy_log1p.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${normalized ? "--normalized ${normalized}" : ''} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.log1p.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
