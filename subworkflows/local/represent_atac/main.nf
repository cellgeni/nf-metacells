nextflow.enable.types = true

include { ATAC_TFIDF } from '../../../modules/cellgeni/atac/tfidf'
include { ATAC_LSI } from '../../../modules/cellgeni/atac/lsi'
include { SNAPATAC2_SELECTFEATURES } from '../../../modules/cellgeni/snapatac2/selectfeatures'
include { SNAPATAC2_SPECTRAL } from '../../../modules/cellgeni/snapatac2/spectral'
include { AssayState; StepOutput } from '../../../workflows/types'

/*
 * An embedding of ATAC assay states, one module per step: TF-IDF -> LSI (muon's procedure, as the v1 modules) for
 * embedding 'lsi', or snapATAC2 feature selection -> spectral embedding for 'spectral'. Counts come from counts_key
 * of each state's object. Emits the embedding output of each state and every step output, for publishing.
 */
workflow REPRESENT_ATAC {
    take:
    states: Channel<AssayState>
    embedding: String
    counts_key: String

    main:
    lsi_states = states.filter { _st -> embedding == 'lsi' }
    tfidf = ATAC_TFIDF(
        lsi_states.map { st -> record(id: st.id, adata: st.adata, counts_key: counts_key, output_prefix: null) }
    )
    lsi = ATAC_LSI(
        lsi_states.join(tfidf, by: 'id')
            .map { st -> record(id: st.id, adata: st.adata, tfidf: st.tfidf, output_prefix: null) }
    )

    spectral_states = states.filter { _st -> embedding == 'spectral' }
    selected = SNAPATAC2_SELECTFEATURES(
        spectral_states.map { st -> record(id: st.id, adata: st.adata, counts_key: counts_key, output_prefix: null) }
    )
    spectral = SNAPATAC2_SPECTRAL(
        spectral_states.join(selected, by: 'id')
            .map { st -> record(id: st.id, adata: st.adata, selected: st.selected, counts_key: counts_key, output_prefix: null) }
    )

    outputs = lsi_states.join(lsi, by: 'id').map { st -> atacOutput(st, 'lsi', st.lsi) }
        .mix(spectral_states.join(spectral, by: 'id').map { st -> atacOutput(st, 'spectral', st.spectral) })
    steps = lsi_states.join(tfidf, by: 'id').map { st -> atacOutput(st, 'tfidf', st.tfidf) }
        .mix(spectral_states.join(selected, by: 'id').map { st -> atacOutput(st, 'selectfeatures', st.selected) })
        .mix(outputs)

    emit:
    embeddings: Channel<StepOutput> = outputs
    steps: Channel<StepOutput> = steps
}

def atacOutput(st: AssayState, step: String, output: Path) -> StepOutput {
    return record(id: "${st.id}.${step}", sample: st.sample, assay: st.assay, step: step, output: output)
}
