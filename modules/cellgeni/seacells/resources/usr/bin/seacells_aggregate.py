#!/usr/bin/env python3

import os
import pickle
import logging
import argparse
from typing import Optional
import muon
import SEACells
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def _axis_sum(matrix, axis: int) -> np.ndarray:
    """Return matrix sums as a 1D ndarray for dense, sparse and matrix inputs."""
    summed = matrix.sum(axis=axis)
    if hasattr(summed, "A1"):
        return summed.A1
    return np.asarray(summed).ravel()


def _safe_n_components(requested: int, n_obs: int, n_vars: int, label: str) -> int:
    """Cap PCA/LSI dimensionality to a value supported by the input shape."""
    max_components = min(n_obs, n_vars) - 1
    if max_components < 1:
        raise ValueError(
            f"Cannot compute {label}: need at least 2 observations and 2 features, "
            f"got n_obs={n_obs}, n_vars={n_vars}"
        )
    effective = min(requested, max_components)
    if effective != requested:
        logging.warning(
            "Reducing %s components from %s to %s for data with n_obs=%s and n_vars=%s",
            label,
            requested,
            effective,
            n_obs,
            n_vars,
        )
    return effective


def _positive_neighbor_counts(data: np.ndarray, n_neighbors: int) -> np.ndarray:
    """Return non-zero distance counts in Scanpy's neighbour graph.

    Palantir's bundled diffusion-map code indexes the 10th non-zero neighbour
    distance from this exact graph. Rows with too few non-zero distances occur
    when many cells have identical low-dimensional coordinates; scipy sparse
    matrices do not store zero-distance neighbour edges.
    """
    temp = sc.AnnData(data)
    sc.pp.neighbors(temp, n_pcs=0, n_neighbors=n_neighbors)
    distances = temp.obsp["distances"].tocsr()
    return np.diff(distances.indptr)


def stabilize_embedding_for_neighbors(
    adata: sc.AnnData,
    components_key: str,
    seacells_n_neighbors: int = 15,
    palantir_knn: int = 30,
    random_seed: int = 0,
) -> sc.AnnData:
    """Break exact embedding ties that make SEACells/Palantir kNN rows empty.

    SEACells and the Palantir version pinned in the image both derive adaptive
    bandwidths from non-zero entries in a sparse kNN distance matrix. If many
    cells have identical PCA/LSI coordinates, their nearest-neighbour distances
    are exactly zero and are dropped by the sparse matrix. Palantir then raises
    `IndexError: index 9 is out of bounds ... size 0`, while SEACells can build
    zero radii and emit divide-by-zero warnings.

    The biological information in exactly tied coordinates is indistinguishable
    to these graph builders, so a tiny deterministic jitter is used only when a
    preflight Scanpy neighbour graph shows too few positive distances.
    """
    if components_key not in adata.obsm:
        raise ValueError(f"{components_key!r} not found in adata.obsm")

    embedding = np.asarray(adata.obsm[components_key])
    if embedding.ndim != 2:
        raise ValueError(
            f"{components_key!r} must be a 2D embedding, got shape {embedding.shape}"
        )
    if not np.all(np.isfinite(embedding)):
        bad_values = int((~np.isfinite(embedding)).sum())
        raise ValueError(
            f"{components_key!r} contains {bad_values} non-finite values; "
            "cannot construct SEACells neighbour graph"
        )

    if adata.n_obs <= 2:
        return adata

    n_neighbors = min(palantir_knn, adata.n_obs - 1)
    if n_neighbors < 2:
        return adata

    # Palantir uses adaptive_k=floor(knn/3), with knn=30 by default, and
    # SEACells uses k//2 for its default k=15. Require enough positive sparse
    # distances for both code paths.
    required_positive = max(
        1,
        min(n_neighbors, palantir_knn // 3),
        min(n_neighbors, seacells_n_neighbors // 2),
    )

    positive_counts = _positive_neighbor_counts(embedding, n_neighbors)
    too_few = positive_counts < required_positive
    if not np.any(too_few):
        return adata

    finite_values = embedding[np.isfinite(embedding)]
    nonzero_std = np.std(embedding, axis=0)
    nonzero_std = nonzero_std[np.isfinite(nonzero_std) & (nonzero_std > 0)]
    if nonzero_std.size:
        scale = float(np.median(nonzero_std))
    elif finite_values.size and np.nanmax(np.abs(finite_values)) > 0:
        scale = float(np.nanmax(np.abs(finite_values)))
    else:
        scale = 1.0

    jitter_scale = max(scale * 1e-6, np.finfo(float).eps)
    rng = np.random.RandomState(random_seed)
    jittered = embedding.astype(np.float64, copy=True)
    jittered += rng.normal(loc=0.0, scale=jitter_scale, size=jittered.shape)
    adata.obsm[components_key] = jittered

    logging.warning(
        "Added deterministic jitter with scale %.3g to %s because %s/%s cells "
        "had fewer than %s positive neighbour distances; this avoids empty "
        "kNN rows in SEACells/Palantir waypoint initialisation",
        jitter_scale,
        components_key,
        int(too_few.sum()),
        adata.n_obs,
        required_positive,
    )
    return adata


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Script validates sample and annotation tables and splits annotation table into separate celltypes"
    )
    parser.add_argument(
        "--adata",
        type=str,
        metavar="<file>",
        help="Specify a path to AnnData object",
    )
    parser.add_argument(
        "--sample",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample name for the file",
    )
    parser.add_argument(
        "--n_metacells",
        metavar="<int>",
        type=int,
        help="Specify a number of metacells to be produced by SEACells",
    )
    parser.add_argument(
        "--gamma",
        type=int,
        metavar="<file>",
        help="Specify a parameter gamma to calculate number of metacells. \
        So n_metacells = n_cells / gamma. Mutually exclusive with n_metacells",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        metavar="<dir>",
        help="Specify an output directory to save results",
    )
    parser.add_argument(
        "--type",
        type=str,
        metavar="<str>",
        choices=["gex", "atac"],
        help="Specify data type",
    )
    parser.add_argument(
        "--n_top_genes",
        type=int,
        metavar="<int>",
        help="Specify number of top genes",
        default=2000,
    )
    parser.add_argument(
        "--n_components",
        type=int,
        metavar="<int>",
        help="Specify number of components to calculate for PCA and SVD",
        default=50,
    )
    parser.add_argument(
        "--celltype_label",
        type=str,
        metavar="<str>",
        help="Specify celltype label",
        default=None,
    )
    parser.add_argument(
        "--convergence_epsilon",
        type=float,
        metavar="<float>",
        default=1e-5,
        help="Specify convergence epsilon",
    )
    parser.add_argument(
        "--min_iter",
        type=int,
        metavar="<int>",
        default=10,
        help="Specify minimum number of iterations",
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        metavar="<int>",
        default=50,
        help="Specify maximum number of iterations",
    )
    parser.add_argument(
        "--n_waypoint_eigs",
        type=int,
        metavar="<int>",
        default=10,
        help="Specify number of components to use for initialization",
    )
    parser.add_argument(
        "--use_sparse",
        action="store_true",
        help="Specify whether to use sparse matrix",
    )
    parser.add_argument(
        "--precomputed",
        type=str,
        metavar="<obsm_key>",
        default=None,
        help="Specify obsm key with precomputed embedding",
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample suffix for barcode if needed",
    )

    return parser


def process_gex(
    adata: sc.AnnData, n_top_genes: int = 2000, n_components: int = 50
) -> sc.AnnData:
    """
    Preprocess GEX data
    Args:
        adata (sc.AnnData): AnnData object
        n_top_genes (int): Number of highly variable genes to select
    Returns:
        sc.AnnData: Preprocessed AnnData object
    """
    logging.info("Starting GEX data preprocessing")
    # normalize and take a logarithm
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    n_highly_variable = int(adata.var["highly_variable"].sum())
    if n_highly_variable < 2:
        raise ValueError(
            f"Cannot compute PCA: only {n_highly_variable} highly variable genes selected"
        )

    # compute PCA
    n_components = _safe_n_components(
        n_components, adata.n_obs, n_highly_variable, "PCA"
    )
    sc.tl.pca(adata, n_comps=n_components, use_highly_variable=True)
    logging.info("Completed GEX data preprocessing")
    return adata


def process_atac(adata: sc.AnnData, n_components: int = 50) -> sc.AnnData:
    """
    Preprocess ATAC data
    Args:
        adata (sc.AnnData): AnnData object
    Returns:
        sc.AnnData: Preprocessed AnnData object
    """
    logging.info("Starting ATAC data preprocessing")

    # muon.atac.pp.tfidf divides by per-feature counts; all-zero peaks create
    # infinite IDF values and noisy RuntimeWarnings, so remove them first.
    feature_counts = _axis_sum(adata.X, axis=0)
    keep_features = np.isfinite(feature_counts) & (feature_counts > 0)
    if not np.all(keep_features):
        logging.warning(
            "Removing %s ATAC features with zero or non-finite total counts before TF-IDF",
            int((~keep_features).sum()),
        )
        adata = adata[:, keep_features].copy()
    if adata.n_vars < 2:
        raise ValueError(
            f"Cannot compute LSI: only {adata.n_vars} non-zero ATAC features remain"
        )

    # compute TF-IDF
    muon.atac.pp.tfidf(adata, scale_factor=1e4)

    # compute LSI
    n_components = _safe_n_components(n_components, adata.n_obs, adata.n_vars, "LSI")
    muon.atac.tl.lsi(adata, n_comps=n_components)
    logging.info("Completed ATAC data preprocessing")
    return adata


def get_metacell_number(
    adata: sc.AnnData, n_metacells: Optional[int], gamma: Optional[int]
) -> int:
    """
    Calculate number of metacells
    Args:
        adata (sc.AnnData): AnnData object
        n_metacells (int): Number of metacells
        gamma (int): Parameter gamma to calculate number of metacells
    Returns:
        int: Number of metacells
    """
    logging.info("Calculating number of metacells")
    if n_metacells is not None and gamma is not None:
        raise ValueError(
            "Both n_metacells and gamma cannot be specified at the same time"
        )
    if n_metacells is None and gamma is None:
        raise ValueError("Either n_metacells or gamma should be specified")
    if gamma is not None and gamma <= 0:
        raise ValueError("gamma must be a positive integer")

    if n_metacells is None:
        n_metacells = max(1, round(adata.n_obs / gamma))

    if n_metacells < 1:
        raise ValueError("n_metacells must be at least 1")
    if n_metacells > adata.n_obs:
        raise ValueError(
            f"n_metacells ({n_metacells}) cannot exceed the number of cells ({adata.n_obs})"
        )
    return int(n_metacells)


def fit_seacells_model(
    adata: sc.AnnData,
    n_metacells: int,
    components_key: str,
    n_waypoint_eigs: int,
    convergence_epsilon: float,
    mit_iter: int,
    max_iter: int,
    use_sparse: bool = False,
) -> SEACells.core.SEACells:
    """
    Fit SEACells model
    Args:
        adata (sc.AnnData): AnnData object
        n_metacells (int): Number of metacells
        components_key (str): Key for components
        n_waypoint_eigs (int): Number of waypoint eigenvectors
        convergence_epsilon (float): Convergence epsilon
        mit_iter (int): Minimum number of iterations
        max_iter (int): Maximum number of iterations
    Returns:
        SEACells.core.SEACells: SEACells model
    """
    logging.info("Fitting SEACells model")

    # Break exact PCA/LSI ties before SEACells constructs its kernel or Palantir
    # computes waypoint diffusion maps. Without this, rows whose neighbours all
    # have zero distance can disappear from sparse distance graphs.
    adata = stabilize_embedding_for_neighbors(adata, components_key)

    # Palantir's waypoint sampler works on diffusion components returned by
    # determine_multiscale_space(). With the Palantir version bundled in this
    # image, n_eigs=1 yields zero diffusion components and crashes with
    # ZeroDivisionError. Conversely, too many diffusion components relative to
    # num_waypoints gives int(num_waypoints / n_components) == 0 and crashes
    # with IndexError. Keep the effective value in the safe range:
    #     2 <= n_waypoint_eigs <= n_metacells + 1
    effective_n_waypoint_eigs = min(n_waypoint_eigs, n_metacells + 1)
    effective_n_waypoint_eigs = max(2, effective_n_waypoint_eigs)
    if effective_n_waypoint_eigs != n_waypoint_eigs:
        logging.warning(
            "Changing n_waypoint_eigs from %s to %s because n_metacells=%s",
            n_waypoint_eigs,
            effective_n_waypoint_eigs,
            n_metacells,
        )

    # create a model
    model = SEACells.core.SEACells(
        adata,
        build_kernel_on=components_key,
        n_SEACells=n_metacells,
        n_waypoint_eigs=effective_n_waypoint_eigs,
        convergence_epsilon=convergence_epsilon,
        use_sparse=use_sparse,
    )

    # construct kernel
    model.construct_kernel_matrix()

    # initialize archetypes
    model.initialize_archetypes()

    # fit model
    model.fit(min_iter=mit_iter, max_iter=max_iter)
    logging.info("Completed fitting SEACells model")
    return model


def plot_assignments(model: SEACells.core.SEACells, output_dir: str):
    """
    Plot assignments
    Args:
        adata (sc.AnnData): AnnData object
        model (sc.AnnData): SEACells model
    """
    logging.info("Plotting assignments")
    # create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 15), gridspec_kw={"wspace": 0.3})

    assignment_weights = model.A_.T
    if hasattr(assignment_weights, "toarray"):
        assignment_weights = assignment_weights.toarray()
    assignment_weights = np.asarray(assignment_weights)

    # non-trivial assignments
    sns.histplot((assignment_weights > 0.1).sum(axis=1), kde=False, ax=ax1)
    ax1.set_title("Non-trivial (> 0.1) assignments per cell")
    ax1.set_xlabel("# Non-trivial SEACell Assignments")
    ax1.set_ylabel("# Cells")

    # weights
    top_n = min(5, assignment_weights.shape[1])
    b = np.partition(assignment_weights, -top_n, axis=1)
    sns.heatmap(np.sort(b[:, -top_n:], axis=1)[:, ::-1], cmap="viridis", vmin=0, ax=ax2)
    ax2.set_title(f"Strength of top {top_n} strongest assignments")
    ax2.set_xlabel("$n^{th}$ strongest assignment")
    fig.savefig(os.path.join(output_dir, "assignments.pdf"))
    plt.close(fig)
    logging.info("Completed plotting assignments")


def get_soft_assignments(model: SEACells.core.SEACells, max_assignments: int = 5):
    """Compute top soft assignments without assuming at least five metacells.

    SEACells.get_soft_assignments() always asks for five assignments. For
    samples with fewer than five metacells it repeats exhausted columns and can
    emit -1 weights. This helper returns min(5, n_metacells) valid assignments.
    """
    assignment_weights = model.A_.T
    if hasattr(assignment_weights, "toarray"):
        assignment_weights = assignment_weights.toarray()
    assignment_weights = np.asarray(assignment_weights)

    if assignment_weights.ndim != 2 or assignment_weights.shape[1] < 1:
        raise ValueError(
            f"Invalid SEACells assignment matrix shape: {assignment_weights.shape}"
        )

    top_n = min(max_assignments, assignment_weights.shape[1])
    order = np.argsort(-assignment_weights, axis=1)[:, :top_n]
    weights = np.take_along_axis(assignment_weights, order, axis=1)

    try:
        archetype_labels = np.asarray(model.get_hard_archetypes())
    except Exception:
        archetype_labels = np.asarray(
            [f"SEACell-{i}" for i in range(assignment_weights.shape[1])]
        )
    labels = archetype_labels[order]

    soft_labels = pd.DataFrame(labels, index=model.ad.obs_names)
    return soft_labels, weights


def compute_celltype_purity(adata: sc.AnnData, celltype_label: str) -> pd.DataFrame:
    """Compute SEACell purity while tolerating missing labels.

    SEACells.evaluate.compute_celltype_purity assumes every SEACell has at least
    one non-missing label. Real annotations can contain all-missing metacells; in
    that case the purity is undefined and is reported as NaN instead of raising.
    """
    rows = []
    for seacell, labels in adata.obs.groupby("SEACell", observed=True)[celltype_label]:
        counts = labels.dropna().value_counts()
        counts = counts[counts > 0]
        if counts.empty:
            rows.append(
                {
                    "SEACell": seacell,
                    celltype_label: pd.NA,
                    f"{celltype_label}_purity": np.nan,
                }
            )
        else:
            rows.append(
                {
                    "SEACell": seacell,
                    celltype_label: counts.index[0],
                    f"{celltype_label}_purity": counts.iloc[0] / counts.sum(),
                }
            )

    return pd.DataFrame(rows).set_index("SEACell")


def _mark_axis_unavailable(ax, title: str, reason: str):
    """Render a placeholder panel for QC metrics that cannot be computed."""
    ax.set_title(title)
    ax.text(0.5, 0.5, reason, ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])


def plot_metacell_stats(
    adata: sc.AnnData,
    output_dir: str,
    components_key: str,
    celltype_label: Optional[str] = None,
):
    """
    Plot metacells
    Args:
        adata (sc.AnnData): AnnData object
        output_dir (str): Output directory
        components_key (str): Key for components
        celltype_label (str): Celltype label
    """
    logging.info("Plotting metacell stats")
    # create plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2, 2, figsize=(10, 15), gridspec_kw={"hspace": 0.3, "wspace": 0.3}
    )

    # metacell sizes
    label_df = adata.obs[["SEACell"]].reset_index()
    sns.histplot(label_df.groupby("SEACell", observed=True).count().iloc[:, 0], ax=ax1)
    ax1.set_title("Metacell Sizes")
    ax1.set_xlabel("# Cells per Metacell")

    # compactness
    try:
        compactness = SEACells.evaluate.compactness(adata, components_key)
        sns.boxplot(data=compactness, y="compactness", ax=ax2)
        ax2.set_title("Compactness")
    except Exception as exc:  # QC should not invalidate completed assignments
        logging.warning("Skipping compactness plot: %s", exc)
        _mark_axis_unavailable(ax2, "Compactness", "unavailable")

    # separation
    try:
        separation = SEACells.evaluate.separation(adata, components_key, nth_nbr=1)
        sns.boxplot(data=separation, y="separation", ax=ax3)
        ax3.set_title("Separation")
    except Exception as exc:  # QC should not invalidate completed assignments
        logging.warning("Skipping separation plot: %s", exc)
        _mark_axis_unavailable(ax3, "Separation", "unavailable")

    # purity
    if celltype_label:
        if celltype_label not in adata.obs.columns:
            logging.warning(
                "Skipping celltype purity plot: %s not found in adata.obs",
                celltype_label,
            )
            _mark_axis_unavailable(ax4, "Celltype Purity", "label not found")
        else:
            purity = compute_celltype_purity(adata, celltype_label)
            purity_col = f"{celltype_label}_purity"
            if purity[purity_col].dropna().empty:
                logging.warning(
                    "Skipping celltype purity plot: no non-missing labels found for %s",
                    celltype_label,
                )
                _mark_axis_unavailable(ax4, "Celltype Purity", "no labels")
            else:
                sns.boxplot(data=purity, y=purity_col, ax=ax4)
                ax4.set_title(f"{celltype_label} Purity")
    else:
        ax4.remove()

    # save plots
    fig.savefig(os.path.join(output_dir, "metacell_stats.pdf"))
    plt.close(fig)
    logging.info("Completed plotting metacell stats")


def evaluate_results(
    adata: sc.AnnData,
    model: SEACells.core.SEACells,
    output_dir: str,
    components_key: str,
    celltype_label: Optional[str] = None,
):
    """
    Evaluate results
    Args:
        adata (sc.AnnData): AnnData object
        model (sc.AnnData): SEACells model
        output_dir (str): Output directory
        components_key (str): Key for components
        celltype_label (str): Celltype label
    """
    logging.info("Evaluating results")
    # plot convergence
    try:
        model.plot_convergence(
            save_as=os.path.join(output_dir, "convergence.pdf"), show=False
        )
        plt.close("all")
    except Exception as exc:  # QC should not invalidate completed assignments
        logging.warning("Skipping convergence plot: %s", exc)

    # plot soft assignments
    try:
        plot_assignments(model, output_dir)
    except Exception as exc:  # QC should not invalidate completed assignments
        logging.warning("Skipping assignments plot: %s", exc)

    # plot metacell stats
    try:
        plot_metacell_stats(
            adata,
            output_dir,
            components_key=components_key,
            celltype_label=celltype_label,
        )
    except Exception as exc:  # QC should not invalidate completed assignments
        logging.warning("Skipping metacell stats plot: %s", exc)
    logging.info("Completed evaluating results")


def format_hard_assignmets(
    assignments: pd.DataFrame, sample: str = None, delimiter: str = None
) -> pd.DataFrame:
    """
    Format hard assignments
    Args:
        assignments (pd.DataFrame): Hard assignments
    Returns:
        pd.DataFrame: Formatted hard assignments
    """
    # change column names
    assignments.columns = ["metacell"]
    assignments.index.name = "barcode"

    # add sample name to index
    if delimiter:
        assignments.index = assignments.index + delimiter + sample
    return assignments


def main():
    """
    Main function
    """
    # parse arguments
    parser = init_parser()
    args = parser.parse_args()

    logging.info("Loading AnnData object")
    # load AnnData object
    adata = sc.read_h5ad(args.adata)
    adata.layers["X_raw"] = adata.X.copy()

    # preprocess data
    logging.info("Preprocessing data of type: %s", args.type)
    if args.precomputed:
        components_key = args.precomputed
        adata_processed = adata.copy()
        if components_key not in adata.obsm.keys():
            raise ValueError("Invalid precomputed key")
    elif args.type == "gex":
        adata_processed = process_gex(adata, args.n_top_genes, args.n_components)
        components_key = "X_pca"
    elif args.type == "atac":
        adata_processed = process_atac(adata, args.n_components)
        components_key = "X_lsi"
    else:
        raise ValueError("Invalid data type")

    # get number of metacells
    n_metacells = get_metacell_number(adata_processed, args.n_metacells, args.gamma)

    # compute metacells
    model = fit_seacells_model(
        adata=adata_processed,
        n_metacells=n_metacells,
        components_key=components_key,
        n_waypoint_eigs=args.n_waypoint_eigs,
        convergence_epsilon=args.convergence_epsilon,
        mit_iter=args.min_iter,
        max_iter=args.max_iter,
        use_sparse=args.use_sparse,
    )

    # assign cells to metacells
    logging.info("Make hard assignments")
    hard_labels = format_hard_assignmets(
        model.get_hard_assignments(), sample=args.sample, delimiter=args.delimiter
    )

    logging.info("Make soft assignments")
    soft_labels, weights = get_soft_assignments(model)

    # save results before optional QC plots so a completed fit is retained even
    # if a plotting/evaluation metric is unavailable for a sample.
    logging.info("Saving results")
    os.makedirs(args.output_dir, exist_ok=True)
    adata_processed.write_h5ad(os.path.join(args.output_dir, "seacell_metacells.h5ad"))
    hard_labels.to_csv(
        os.path.join(args.output_dir, "seacell_metacell_assignments.csv")
    )
    soft_labels.to_csv(os.path.join(args.output_dir, "seacell_soft_assignments.csv"))
    np.save(os.path.join(args.output_dir, "seacell_weights.npy"), weights)
    with open(os.path.join(args.output_dir, "seacell_model.pkl"), "wb") as file:
        pickle.dump(model, file)
    logging.info("Successfully saved results")

    # evaluate results
    evaluate_results(
        adata_processed,
        model,
        args.output_dir,
        components_key=components_key,
        celltype_label=args.celltype_label,
    )


if __name__ == "__main__":
    main()
