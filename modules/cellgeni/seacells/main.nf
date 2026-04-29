/*
 * Module: cellgeni/seacells
 */

process SEACELLS {
    tag "${meta.id}"
    container 'docker://quay.io/cellgeni/seacells:latest'

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.csv"), emit: csv
    tuple val(meta), path("*.npy"), emit: npy
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.pkl"), emit: pkl
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: "--type gex --n_top_genes 2000 --n_components 50 --convergence_epsilon 0.00001 --min_iter 10 --max_iter 50"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seacells_aggregate.py \
        ${args} \
        --adata ${adata} \
        --sample ${prefix} \
        --output_dir .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | awk '{print \$2}' )
        SEACells: \$( python -c "import SEACells; print(SEACells.__version__)" )
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
        muon: \$( python -c "import muon; print(muon.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.csv
    touch ${prefix}.npy
    touch ${prefix}.pdf
    touch ${prefix}.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | awk '{print \$2}' )
        SEACells: \$( python -c "import SEACells; print(SEACells.__version__)" )
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
        muon: \$( python -c "import muon; print(muon.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """
}
