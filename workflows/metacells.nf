nextflow.enable.types = true

include { INGEST } from '../subworkflows/local/ingest'
include { REPRESENT_RNA } from '../subworkflows/local/represent_rna'
include { REPRESENT_ATAC } from '../subworkflows/local/represent_atac'
include { SEACELLS } from '../modules/cellgeni/seacells'
include { HMETACELLS } from '../modules/cellgeni/hmetacells'
include { Sample; AssayState; StepOutput; Embedding; Membership; MethodFiles } from './types'

/*
 * Ingest every assay of each sample, compute the embedding of the sample's membership assay (the
 * --membership_assay when the sample has it, otherwise its only assay) unless --embedding_key names one in the
 * object, then run each selected method on that embedding: PCA for RNA, LSI or spectral (--atac_embedding) for ATAC.
 *
 * Both methods are 2.0 modules: the embedding step output goes in and a `barcode,metacell` table comes out,
 * beside a JSON summary. They see only the membership assay; aggregation over the sample's other assays is the
 * next step of phase 5 (docs/architecture.md).
 */
workflow METACELLS {
    take:
    samples: Channel<Sample>
    methods: List<String>
    membership_assay: String
    cell_metadata: Path?
    barcode_delimiter: String?
    embedding_key: String?
    atac_embedding: String
    rna_counts_key: String
    atac_counts_key: String

    main:
    ingest = INGEST(samples, cell_metadata, barcode_delimiter)

    bases = samples.map { s -> record(sample: s.id, basis: basisAssay(s, membership_assay)) }
    basis = ingest.states.join(bases, by: 'sample').filter { st -> st.assay == st.basis }

    compute = basis.filter { _st -> !embedding_key }
    rna = REPRESENT_RNA(compute.filter { st -> st.assay == 'rna' }, rna_counts_key)
    atac = REPRESENT_ATAC(compute.filter { st -> st.assay == 'atac' }, atac_embedding, atac_counts_key)
    computed = compute.join(
        rna.embeddings.mix(atac.embeddings)
            .map { o -> record(id: "${o.sample}.${o.assay}", step: o.step, output: o.output) },
        by: 'id'
    )
    // Both methods read the embedding straight from the step output: no assembled copy stands between them
    embeddings = computed
        .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: "obsm/X_${st.step}",
                            adata: st.adata, embedding: st.output) }
        .mix(
            basis.filter { _st -> embedding_key }
                .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: "obsm/${embedding_key}",
                                    adata: st.adata, embedding: null) }
        )

    seacells_in = embeddings.filter { _e -> 'seacells' in methods }
    seacells = seacells_in.join(
        SEACELLS(
            seacells_in.map { e -> record(id: e.id, adata: e.adata, embedding: e.embedding,
                                          embedding_key: e.key, output_prefix: null) }
        ),
        by: 'id'
    )

    hmetacells_in = embeddings.filter { _e -> 'hmetacells' in methods }
    hmetacells = hmetacells_in.join(
        HMETACELLS(
            hmetacells_in.map { e -> record(id: e.id, adata: e.adata, embedding: e.embedding,
                                            embedding_key: e.key, output_prefix: null) }
        ),
        by: 'id'
    )

    memberships = seacells
        .map { r -> record(id: "${r.id}.seacells", sample: r.sample, assay: r.assay, method: 'seacells',
                           membership: r.membership) }
        .mix(
            hmetacells.map { r -> record(id: "${r.id}.hmetacells", sample: r.sample, assay: r.assay,
                                         method: 'hmetacells', membership: r.membership) }
        )

    extras = seacells
        .map { r -> record(id: "${r.id}.seacells.summary", sample: r.sample, method: 'seacells',
                           files: [r.summary]) }
        .mix(
            hmetacells.map { r -> record(id: "${r.id}.hmetacells.summary", sample: r.sample, method: 'hmetacells',
                                         files: [r.summary]) }
        )

    emit:
    states: Channel<AssayState> = ingest.states
    written: Channel<AssayState> = ingest.written
    representations: Channel<StepOutput> = rna.steps.mix(atac.steps)
    embeddings: Channel<Embedding> = embeddings
    memberships: Channel<Membership> = memberships
    extras: Channel<MethodFiles> = extras
}

def basisAssay(s: Sample, preferred: String) -> String {
    if (preferred == 'atac') {
        return s.atac ? 'atac' : 'rna'
    }
    return s.rna ? 'rna' : 'atac'
}

