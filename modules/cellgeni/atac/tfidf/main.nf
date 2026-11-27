nextflow.enable.types = true

/*
 * TF-IDF of peak counts into layers/tfidf, with muon's formula (the v1 modules' ATAC preprocessing): NumPy/SciPy
 * on CPU, CuPy when the task has an accelerator. ext.args: --scale-factor <float> (default 1e4), --no-log-tf,
 * --no-log-idf, --log-tfidf (with both previous). Peaks and cells without counts stay zero; no peak is removed.
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process ATAC_TFIDF {
    tag id
    label 'process_gpu'
    container task.accelerator ? 'quay.io/cellgeni/rapids-singlecell:0.18.1' : 'quay.io/cellgeni/scanpy:1.12.6'

    input:
    record(id: String, adata: Path, counts_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, tfidf: file("${output_prefix ?: id}.tfidf.*"))

    topic:
    record(
        process: task.process,
        tool: task.accelerator ? 'cupy' : 'scipy',
        version: eval("scio.py version ${task.accelerator ? 'cupy-cuda12x' : 'scipy'}")
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    atac_tfidf.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${counts_key ? "--counts-key ${counts_key}" : ''} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.tfidf.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
