/*
 * Module: cellgeni/hmetacells
 */

process HMETACELLS {
    tag "${meta.id}"
    container 'quay.io/cellgeni/metacells-python:latest'

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: "--type gex --n_min 5 --n_max 20 --method ward --n_top_genes 2000 --n_components 50 --n_neighbors 15"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hierarchial_metacells.py \
        ${args} \
        --adata ${adata} \
        --sample ${prefix} \
        --output hierarchial_metacells.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | awk '{print \$2}' )
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
        muon: \$( python -c "import muon; print(muon.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    touch hierarchial_metacells.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version 2>&1 | awk '{print \$2}' )
        anndata: \$( python -c "import anndata; print(anndata.__version__)" )
        scanpy: \$( python -c "import scanpy; print(scanpy.__version__)" )
        muon: \$( python -c "import muon; print(muon.__version__)" )
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """
}