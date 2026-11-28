nextflow.enable.types = true

include { INGEST } from '../subworkflows/local/ingest'
include { REPRESENT_RNA } from '../subworkflows/local/represent_rna'
include { REPRESENT_ATAC } from '../subworkflows/local/represent_atac'
include { ANNDATA_ASSEMBLE } from '../modules/cellgeni/anndata/assemble'
include { SEACELLS } from '../modules/cellgeni/seacells'
include { HMETACELLS } from '../modules/cellgeni/hmetacells'
include { Sample; AssayState; StepOutput; Embedding; Membership; MethodFiles } from './types'

/*
 * Ingest every assay of each sample, compute the embedding of the sample's membership assay (the
 * --membership_assay when the sample has it, otherwise its only assay) unless --embedding_key names one in the
 * object, then run each selected method on that embedding: PCA for RNA, LSI or spectral (--atac_embedding) for ATAC.
 *
 * SEACells is the 2.0 module: the embedding step output goes in, a `barcode,metacell` table comes out. hmetacells
 * is still the v1 module until the rest of phase 5 (docs/architecture.md): it reads the embedding from obsm of the
 * object it aggregates (--precomputed), so the step output is assembled with the base object first. Both see only
 * the membership assay.
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
    assembled = ANNDATA_ASSEMBLE(
        computed.map { st -> record(id: st.id, adata: st.adata, sources: [st.output], output_prefix: null) }
    )
    embeddings = computed.join(assembled, by: 'id')
        .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: "X_${st.step}", adata: st.assembled) }
        .mix(
            basis.filter { _st -> embedding_key }
                .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: embedding_key, adata: st.adata) }
        )
    // SEACells 2.0 reads the embedding straight from the step output, so it does not wait for the assembled copy
    method_input = computed
        .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: "X_${st.step}",
                            adata: st.adata, embedding: st.output) }
        .mix(
            basis.filter { _st -> embedding_key }
                .map { st -> record(id: st.id, sample: st.sample, assay: st.assay, key: embedding_key,
                                    adata: st.adata, embedding: null) }
        )

    seacells_in = method_input.filter { _st -> 'seacells' in methods }
    seacells = seacells_in.join(
        SEACELLS(
            seacells_in.map { e -> record(id: e.id, adata: e.adata, embedding: e.embedding,
                                          embedding_key: e.key, output_prefix: null) }
        ),
        by: 'id'
    )

    // v1 hmetacells still reads the embedding from obsm of the assembled object (--precomputed)
    v1_input = embeddings.map { e -> tuple([id: e.sample, assay: e.assay, embedding_key: e.key], e.adata) }
    hmetacells = HMETACELLS(v1_input.filter { _st -> 'hmetacells' in methods })

    memberships = seacells
        .map { r -> record(id: "${r.id}.seacells", sample: r.sample, assay: r.assay, method: 'seacells',
                           membership: r.membership) }
        .mix(hmetacells.csv.map { meta, table -> membership(meta, 'hmetacells', table) })

    extras = seacells.map { r -> record(id: "${r.id}.seacells.summary", sample: r.sample, method: 'seacells',
                                        files: [r.summary]) }

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

def membership(meta: Map, method: String, table: Path) -> Membership {
    return record(id: "${meta.id}.${meta.assay}.${method}", sample: meta.id, assay: meta.assay, method: method, membership: table)
}
