process HierarchialAggregate {
    tag "Running SEAcells metacells aggregation foor ${item}"
    input:
        tuple val(item), path(adata)
        val(n_min)
        val(n_max)
        val(type)
        val(n_top_genes)
        val(n_components)
        val(celltype_label)
        val(n_neighbors)
        val(precomputed)
        val(method)
        val(sample_prefix)
        val(prefix_delimiter)
    output:
        path("${item}")
    script:
        """
        mkdir -p ${item}
        hierarchial_metacells.py \
            --adata ${adata} \
            --celltype_label ${celltype_label} \
            --method ${method} \
            --output metacells.csv \
            --n_min ${n_min} \
            --n_max ${n_max} \
            --type ${type} \
            --n_top_genes ${n_top_genes} \
            --n_components ${n_components} \
            --n_neighbors ${n_neighbors} \
            ${precomputed ? "--precomputed ${precomputed}" : ""} \
            ${sample_prefix ? "--sample_prefix ${sample_prefix}" : ""} \
            ${prefix_delimiter ? "--prefix_delimiter ${prefix_delimiter}" : ""}
        """
}