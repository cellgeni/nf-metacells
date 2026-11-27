nextflow.enable.types = true

/*
 * Marks informative peaks or tiles in var/selected with snapATAC2 select_features; never subsets: snapATAC2 on
 * CPU. ext.args: --n-features <int> (default 500000), --filter-lower-quantile, --filter-upper-quantile (default
 * 0.005), --whitelist <bed>, --blacklist <bed> (need region var names), --max-iter <int> (default 1).
 * Output (ext.output_mode): raw directory (default), slim or attach object (ext.adata_format: h5ad or zarr);
 * sparse matrices in raw directories as ext.sparse_format: npz (default) or mtx.
 */
process SNAPATAC2_SELECTFEATURES {
    tag id
    label 'process_single'
    container 'quay.io/cellgeni/snapatac2:2.10.0'

    input:
    record(id: String, adata: Path, counts_key: String?, output_prefix: String?)

    stage:
    stageAs adata, 'input/*'

    output:
    record(id: id, selected: file("${output_prefix ?: id}.selectfeatures.*"))

    topic:
    record(
        process: task.process,
        tool: 'snapatac2',
        version: eval('scio.py version snapatac2')
    ) >> 'versions'

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} NUMBA_NUM_THREADS=${task.cpus}
    snapatac2_selectfeatures.py \\
        --adata ${adata} \\
        --prefix ${output_prefix ?: id} \\
        ${counts_key ? "--counts-key ${counts_key}" : ''} \\
        --threads ${task.cpus} \\
        --output-mode ${task.ext.output_mode ?: 'raw'} \\
        --adata-format ${task.ext.adata_format ?: 'h5ad'} \\
        --sparse-format ${task.ext.sparse_format ?: 'npz'} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def mode = task.ext.output_mode ?: 'raw'
    def out = "${output_prefix ?: id}.selectfeatures.${mode == 'raw' ? 'raw' : (task.ext.adata_format ?: 'h5ad')}"
    """
    ${out.endsWith('.h5ad') ? "touch ${out}" : "mkdir ${out}"}
    """
}
