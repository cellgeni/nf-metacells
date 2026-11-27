// Requests the Seqera Containers (Wave) builds of the conda images and prints their names (see docker/README.md):
//   cd docker/wave && nextflow inspect main.nf -concretize
process SCANPY_IMAGE {
    conda "${projectDir}/../scanpy/environment.yml"

    script:
    'true'
}

process SNAPATAC2_IMAGE {
    conda "${projectDir}/../snapatac2/environment.yml"

    script:
    'true'
}

workflow {
    SCANPY_IMAGE()
    SNAPATAC2_IMAGE()
}
