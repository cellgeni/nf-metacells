process TO_H5AD {
    tag "Converting ${sample_id}'s file to .h5ad"
    input:
        tuple val(sample_id), path(input)
        val(delimiter)
    output:
        tuple val(sample_id), path("${sample_id}.h5ad")
    script:
        """
        convert_to_h5ad.py \
            --input ${input} \
            --sample_id ${sample_id} \
            --delimiter ${delimiter} \
            --output ${sample_id}.h5ad
        """
}

process ATTACH_CELL_METADATA {
    tag "Attaching cell metadata to .obs section of AnnData object for ${sample}"
    input:
        tuple val(sample), path(h5ad)
        val(metadata)
        val(barcode_column)
    output:
        tuple val(sample), path("${sample}.h5ad")
    script:
        """
        attach_annotation.py \
            --h5ad_file ${h5ad} \
            --sample_id ${sample} \
            --metadata ${metadata} \
            --barcode_column ${barcode_column} \
            --output ${sample}.h5ad
        """
}