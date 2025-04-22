#!/usr/bin/env python3

import logging
import argparse
import csv
from typing import Dict
import numpy as np
import muon
import scanpy as sc
from sknetwork.hierarchy import Paris, LouvainHierarchy, cut_balanced

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
        required=True,
    )
    parser.add_argument(
        "--sample",
        type=str,
        metavar="<str>",
        default=None,
        help="Specify sample name for the file",
    )
    parser.add_argument(
        "--n_min",
        metavar="<int>",
        type=int,
        help="Specify a min number of cells in metacell",
        default=10,
    )
    parser.add_argument(
        "--n_max",
        metavar="<int>",
        type=int,
        help="Specify a max number of cells in metacell",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        metavar="<file>",
        help="Specify an output filename to save results",
    )
    parser.add_argument(
        "--type",
        type=str,
        metavar="<gex|atac>",
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
        "--n_neighbors",
        type=int,
        metavar="<int>",
        help="Specify number of nearest neighbors",
        default=15,
    )
    parser.add_argument(
        "--celltype_label",
        type=str,
        metavar="<str>",
        help="Specify celltype label",
        default=None,
    )
    parser.add_argument(
        "--precomputed",
        type=str,
        metavar="<obsm_key>",
        default=None,
        help="Specify obsm key with precomputed embedding",
    )
    parser.add_argument(
        "--method",
        type=str,
        metavar="<paris|louvain>",
        default="louvain",
        help="Specify method to use for cell clustering",
        required=True,
    )
    parser.add_argument(
        "--sample_suffix",
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


def cluster_cells(
    adata: sc.AnnData,
    method: str = "louvain",
    n_components: int = 50,
    n_neighbors: int = 15,
    components_key: str = "X_pca",
) -> np.ndarray:
    """
    Cluster cells
    Args:
        adata (sc.AnnData): AnnData object
        method (str): Method to use for clustering
        n_neighbors (int): Number of nearest neighbors
    Returns:
        np.ndarray: Cluster assignments
    """
    # Calculate nearest neighbors
    logging.info("Calculating nearest neighbors")
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_components,
        use_rep=components_key,
        random_state=4,
    )

    logging.info("Starting cell clustering")
    if method != "paris" and method != "louvain":
        raise ValueError(f"Invalid method: {method}")

    model = LouvainHierarchy() if method == "louvain" else Paris()
    dendrogram = model.fit_transform(adata.obsp["connectivities"])
    logging.info("Completed cell clustering")
    return dendrogram


def calculate_metacells(
    adata: sc.AnnData,
    n_min: int,
    n_max: int,
    celltype_label: str,
    method: str = "louvain",
    n_components: int = 50,
    n_neighbors: int = 15,
    components_key: str = "X_pca",
) -> Dict[str, str]:
    """
    Calculate metacells
    Args:
        adata (sc.AnnData): AnnData object
        n_min (int): Minimum number of cells in metacell
        n_max (int): Maximum number of cells in metacell
        method (str): Method to use for clustering
        n_neighbors (int): Number of nearest neighbors
    Returns:
        Dict[str, str]: Metacell assignments
    """
    metacell_dict = dict()
    for celltype in adata.obs[celltype_label].unique():
        # subset adata by celltype
        logging.info("Calculating metacells for celltype: %s", celltype)
        adata_celltype = adata[adata.obs[celltype_label] == celltype].copy()

        # check if number of cells greater than cluster size
        if adata_celltype.n_obs < n_max:
            adata_celltype.obs["metacell"] = f"{celltype}_0"
            metacell_dict.update(adata_celltype.obs["metacell"].to_dict())
            continue

        # cluster cells
        dendrogram = cluster_cells(
            adata_celltype,
            method=method,
            n_components=n_components,
            n_neighbors=n_neighbors,
            components_key=components_key,
        )

        # cut dendrogram
        clusters = cut_balanced(dendrogram, max_cluster_size=n_max)
        cluster_count = np.bincount(clusters)
        small_clusters = np.where((cluster_count < n_min) & (cluster_count > 0))[0]

        # warn if small clusters
        if small_clusters.size > 0:
            logging.warning("There are small clusters: %s", small_clusters)

        # save metacells to dict
        adata_celltype.obs["metacell"] = f"{celltype}_" + clusters.astype(str)
        metacell_dict.update(adata_celltype.obs["metacell"].to_dict())
    return metacell_dict


def main():
    # parse arguments
    parser = init_parser()
    args = parser.parse_args()

    # load AnnData object
    logging.info("Loading AnnData object")
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

    # calculate metacells
    metacell_dict = calculate_metacells(
        adata_processed,
        args.n_min,
        args.n_max,
        args.celltype_label,
        args.method,
        args.n_components,
        args.n_neighbors,
        components_key,
    )

    # save metacells to csv file
    logging.info("Saving metacells to csv file")
    with open(args.output, mode="w", encoding="utf-8") as csv_file:
        # create a writer object
        writer = csv.DictWriter(csv_file, fieldnames=["barcode", "metacell"])

        # write the data
        writer.writeheader()
        for barcode, metacell in metacell_dict.items():
            barcode = (
                barcode + args.sample_suffix + args.sample
                if args.sample_suffix
                else barcode
            )
            writer.writerow({"barcode": barcode, "metacell": metacell})


if __name__ == "__main__":
    main()
