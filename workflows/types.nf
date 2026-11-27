nextflow.enable.types = true

/*
 * Pipeline records (docs/architecture.md). Every channel item has a unique `id`; items of different kinds are
 * joined on `sample`. Records only grow: steps added in later phases add nullable fields.
 */

// One samplesheet row
record Sample {
    id: String
    rna: Path?
    atac: Path?
    metadata: Path?
}

// One assay of a sample, ready for the methods; id = <sample>.<assay>
record AssayState {
    id: String
    sample: String
    assay: String
    adata: Path
}

// One representation step's output for one assay; id = <sample>.<assay>.<step>
record StepOutput {
    id: String
    sample: String
    assay: String
    step: String
    output: Path
}

// The embedding the methods of a sample run on; id = <sample>.<assay>. `adata` holds obsm[key] and the counts.
record Embedding {
    id: String
    sample: String
    assay: String
    key: String
    adata: Path
}

// One method run on one sample; id = <sample>.<assay>.<method>. `assay` is the assay that drove membership.
record Membership {
    id: String
    sample: String
    assay: String
    method: String
    membership: Path
}

// Files a method writes besides its membership table, published next to it
record MethodFiles {
    id: String
    sample: String
    method: String
    files: List<Path>
}
