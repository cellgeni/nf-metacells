nextflow.enable.types = true

include { SCANPY_NORMALIZETOTAL } from '../../../modules/cellgeni/scanpy/normalizetotal'
include { SCANPY_LOG1P } from '../../../modules/cellgeni/scanpy/log1p'
include { SCANPY_HVG } from '../../../modules/cellgeni/scanpy/hvg'
include { SCANPY_PCA } from '../../../modules/cellgeni/scanpy/pca'
include { AssayState; StepOutput } from '../../../workflows/types'

/*
 * PCA of RNA assay states, one module per step: normalize_total -> log1p -> highly variable genes -> PCA on those
 * genes (the v1 modules' preprocessing). Counts come from counts_key of each state's object. Emits the PCA output
 * of each state (obsm/X_pca) and every step output, for publishing.
 */
workflow REPRESENT_RNA {
    take:
    states: Channel<AssayState>
    counts_key: String

    main:
    normalized = SCANPY_NORMALIZETOTAL(
        states.map { st -> record(id: st.id, adata: st.adata, counts_key: counts_key, output_prefix: null) }
    )
    log1p = SCANPY_LOG1P(
        states.join(normalized, by: 'id')
            .map { st -> record(id: st.id, adata: st.adata, normalized: st.normalized, output_prefix: null) }
    )
    hvg = SCANPY_HVG(
        states.join(log1p, by: 'id')
            .map { st -> record(id: st.id, adata: st.adata, log1p: st.log1p, counts_key: counts_key, output_prefix: null) }
    )
    pca = SCANPY_PCA(
        states.join(log1p, by: 'id').join(hvg, by: 'id')
            .map { st -> record(id: st.id, adata: st.adata, matrix: st.log1p, mask: st.hvg, output_prefix: null) }
    )

    outputs = states.join(pca, by: 'id').map { st -> stepOutput(st, 'pca', st.pca) }
    steps = states.join(normalized, by: 'id').map { st -> stepOutput(st, 'normalizetotal', st.normalized) }
        .mix(states.join(log1p, by: 'id').map { st -> stepOutput(st, 'log1p', st.log1p) })
        .mix(states.join(hvg, by: 'id').map { st -> stepOutput(st, 'hvg', st.hvg) })
        .mix(outputs)

    emit:
    embeddings: Channel<StepOutput> = outputs
    steps: Channel<StepOutput> = steps
}

def stepOutput(st: AssayState, step: String, output: Path) -> StepOutput {
    return record(id: "${st.id}.${step}", sample: st.sample, assay: st.assay, step: step, output: output)
}
