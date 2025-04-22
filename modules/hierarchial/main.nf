process HierarchialAggregate {
    tag "Running hierarchial clustering for ${sample}"
    input:
        tuple val(sample), path(adata)
        val(n_min)
        val(n_max)
        val(type)
        val(n_top_genes)
        val(n_components)
        val(celltype_label)
        val(n_neighbors)
        val(precomputed)
        val(method)
        val(sample_suffix)
    output:
        path("hierarchial_metacells.csv")
    script:
        """
        hierarchial_metacells.py \
            --adata ${adata} \
            --sample ${sample} \
            --celltype_label ${celltype_label} \
            --method ${method} \
            --output hierarchial_metacells.csv \
            --n_min ${n_min} \
            --n_max ${n_max} \
            --type ${type} \
            --n_top_genes ${n_top_genes} \
            --n_components ${n_components} \
            --n_neighbors ${n_neighbors} \
            ${precomputed ? "--precomputed ${precomputed}" : ""} \
            ${sample_suffix ? "--sample_suffix ${sample_suffix}" : ""}
        """
}