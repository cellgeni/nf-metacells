process H5_TO_H5AD {
    tag "Running SEAcells metacells aggregation foor ${sample}"
    input:
        tuple val(sample), path(h5file)
        val(celltype_annotation)
        val(celltype_label)
        val(delimiter)
    output:
        tuple val(sample), path("${sample}.h5ad")
    script:
        """
        h5_to_h5ad.py \
            --h5file ${h5file} \
            --sample ${sample} \
            ${celltype_annotation ? "--celltype_annotation ${celltype_annotation}" : ""} \
            ${celltype_label ? "--celltype_label ${celltype_label}" : ""} \
            ${delimiter ? "--delimiter ${delimiter}" : ""} \
            --output ${sample}.h5ad
        """
}