nextflow.enable.types = true

/*
 * Spectral embedding (Laplacian eigenmaps of cosine similarity) of peak counts into obsm/X_spectral and
 * uns/spectral_eigenvalue: snapATAC2 on CPU, a CuPy port of its matrix-free cosine algorithm when the task has an
 * accelerator. ext.args: --n-comps <int> (default 30; fewer columns when non-positive eigenvalues are dropped),
 * --features-key (default var/selected), --all-features, --distance-metric cosine|jaccard, --sample-size
 * (Nystrom, CPU), --no-weighted-by-sd, --seed (default 0). The GPU path computes the full cosine embedding only.
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SNAPATAC2_SPECTRAL {
    tag id
    label 'process_gpu'
    container task.accelerator ? 'quay.io/cellgeni/rapids-singlecell:0.18.1' : 'quay.io/cellgeni/snapatac2:2.10.0'

    input:
    record(id: String, adata: Path, selected: Path?, counts_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, spectral: file("${output_prefix ?: id}.spectral.*"))

    topic:
    record(
        process: task.process,
        tool: task.accelerator ? 'cupy' : 'snapatac2',
        version: eval("scio.py version ${task.accelerator ? 'cupy-cuda12x' : 'snapatac2'}")
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    snapatac2_spectral.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${selected ? "--selected ${selected}" : ''} \\
        ${counts_key ? "--counts-key ${counts_key}" : ''} \\
        --threads ${task.cpus} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.accelerator ? '--gpu' : ''} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.spectral.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
