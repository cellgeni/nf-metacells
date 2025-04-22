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

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


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

    # compute PCA
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
    # compute TF-IDF
    muon.atac.pp.tfidf(adata, scale_factor=1e4)

    # compute LSI
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
    if n_metacells is not None:
        return n_metacells
    else:
        return round(adata.n_obs / gamma)


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
    # create a model
    model = SEACells.core.SEACells(
        adata,
        build_kernel_on=components_key,
        n_SEACells=n_metacells,
        n_waypoint_eigs=n_waypoint_eigs,
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

    # non-trivial assignments
    sns.displot((model.A_.T > 0.1).sum(axis=1), kde=False, ax=ax1)
    ax1.set_title("Non-trivial (> 0.1) assignments per cell")
    ax1.set_xlabel("# Non-trivial SEACell Assignments")
    ax1.set_ylabel("# Cells")

    # weights
    b = np.partition(model.A_.T, -5)
    sns.heatmap(np.sort(b[:, -5:])[:, ::-1], cmap="viridis", vmin=0, ax=ax2)
    ax2.set_title("Strength of top 5 strongest assignments")
    ax2.set_xlabel("$n^{th}$ strongest assignment")
    plt.savefig(os.path.join(output_dir, "assignments.pdf"))
    logging.info("Completed plotting assignments")


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
    sns.histplot(label_df.groupby("SEACell").count().iloc[:, 0], ax=ax1)
    ax1.set_title("Metacell Sizes")
    ax1.set_xlabel("# Cells per Metacell")

    # compactness
    compactness = SEACells.evaluate.compactness(adata, components_key)
    sns.boxplot(data=compactness, y="compactness", ax=ax2)
    ax2.set_title("Compactness")

    # separation
    separation = SEACells.evaluate.separation(adata, components_key, nth_nbr=1)
    sns.boxplot(data=separation, y="separation", ax=ax3)
    ax3.set_title("Separation")

    # purity
    if celltype_label:
        purity = SEACells.evaluate.compute_celltype_purity(adata, "celltype")
        sns.boxplot(data=purity, y="celltype_purity", ax=ax4)
        ax4.set_title("Celltype Purity")
    else:
        ax4.remove()

    # save plots
    plt.savefig(os.path.join(output_dir, "metacell_stats.pdf"))
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
    model.plot_convergence(
        save_as=os.path.join(output_dir, "convergence.pdf"), show=False
    )

    # plot soft assignments
    plot_assignments(model, output_dir)

    # plot metacell stats
    plot_metacell_stats(
        adata, output_dir, components_key=components_key, celltype_label=celltype_label
    )
    logging.info("Completed evaluating results")


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
    hard_labels = model.get_hard_assignments()
    logging.info("Make soft assignments")
    soft_labels, weights = model.get_soft_assignments()

    # evaluate results
    evaluate_results(
        adata_processed,
        model,
        args.output_dir,
        components_key=components_key,
        celltype_label=args.celltype_label,
    )

    # save results
    logging.info("Saving results")
    adata_processed.write_h5ad(os.path.join(args.output_dir, "seacell_metacells.h5ad"))
    hard_labels.to_csv(
        os.path.join(args.output_dir, "seacell_metacell_assignments.csv")
    )
    soft_labels.to_csv(os.path.join(args.output_dir, "seacell_soft_assignments.csv"))
    np.save(os.path.join(args.output_dir, "seacell_weights.npy"), weights)
    with open(os.path.join(args.output_dir, "seacell_model.pkl"), "wb") as file:
        pickle.dump(model, file)
    logging.info("Successfully saved results")


if __name__ == "__main__":
    main()
