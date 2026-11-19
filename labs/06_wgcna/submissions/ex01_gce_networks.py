"""
Exercise — Gene Co-Expression Networks (GCEs) — Building the network and detecting modules

Objective:
- Build a co-expression network from an RNA-Seq expression matrix
- Detect gene modules (communities) using a Louvain-like algorithm (or alternative)

Instructions (in the lab):
1) Data preparation
   - Download and prepare the expression matrix (e.g., GSE115469) into a CSV with:
     * rows = genes (index), columns = samples (sample IDs)
   - Save the file to: data/work/<handle>/lab06/expression_matrix.csv

2) Preprocessing
   - log2(x + 1)
   - filter low-variance genes

3) Correlation → Adjacency
   - complete the `correlation_matrix` function
   - the function `adjacency_from_correlation` is already implemented

4) Graph + Modules
   - build the graph with NetworkX
   - detect modules (Louvain or alternative)
   - export the gene → module mapping to submissions/<handle>/modules_<handle>.csv

Note:
- Document in <github_handle>_notes.md: correlation metric, threshold, short observations.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------
# Config — fill in as needed
# --------------------------
HANDLE = "AlexZvinca"

INPUT_CSV = Path(
    f"data/work/{HANDLE}/lab06/expression_matrix.csv"
)

OUTPUT_DIR = Path(
    f"labs/06_wgcna/submissions/{HANDLE}"
)

OUTPUT_CSV = OUTPUT_DIR / f"modules_{HANDLE}.csv"

CORR_METHOD = "spearman"   # TODO: "pearson" or "spearman"
VARIANCE_THRESHOLD = 0.5      # threshold for filtering genes
ADJ_THRESHOLD = 0.7           # threshold for |cor| (e.g., 0.6)
USE_ABS_CORR = True            # True => use |cor| when thresholding
MAKE_UNDIRECTED = True         # co-expression networks are usually undirected

# True if the input contains raw, nonnegative RNA-Seq counts.
# False if the data is already library-size normalized.
NORMALIZE_COUNTS = False

# Limit the number of genes to keep the network manageable.
MAX_GENES = 1000

def read_expression_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Expression matrix not found: {path}"
        )

    # First column contains gene identifiers.
    df = pd.read_csv(path, index_col=0)

    if df.empty:
        raise ValueError("Expression matrix is empty.")

    if df.index.has_duplicates:
        raise ValueError("Duplicate gene IDs found.")

    # Convert expression columns to numeric values.
    df = df.apply(pd.to_numeric, errors="raise")

    if df.isna().any().any():
        raise ValueError("Expression matrix contains missing values.")

    if not np.isfinite(df.to_numpy()).all():
        raise ValueError("Expression matrix contains non-finite values.")

    if (df < 0).any().any():
        raise ValueError("Expression values must be nonnegative.")

    if df.shape[1] < 3:
        raise ValueError("At least 3 samples are recommended.")

    return df

def normalize_counts(df: pd.DataFrame) -> pd.DataFrame:
    # Counts per million (CPM) normalization
    library_sizes = df.sum(axis=0)

    if (library_sizes == 0).any():
        raise ValueError("A sample has zero total counts.")

    return df.div(library_sizes, axis=1) * 1_000_000

def log_and_filter(
    df: pd.DataFrame,
    variance_threshold: float
) -> pd.DataFrame:

    # Apply log2(x + 1)
    df_log = np.log2(df + 1)

    # Calculate variance across samples for each gene.
    variances = df_log.var(axis=1)

    # Keep genes with enough variation.
    df_filtered = df_log.loc[
        variances > variance_threshold
    ]

    # Keep the most variable genes if there are too many.
    if len(df_filtered) > MAX_GENES:
        top_genes = (
            variances.loc[df_filtered.index]
            .nlargest(MAX_GENES)
            .index
        )
        df_filtered = df_filtered.loc[top_genes]

    if len(df_filtered) < 2:
        raise ValueError(
            "Too few genes remain. Check the variance threshold."
        )

    return df_filtered


def correlation_matrix(
    df: pd.DataFrame,
    method: str = "spearman",
    use_abs: bool = True
) -> pd.DataFrame:

    if method not in ("pearson", "spearman"):
        raise ValueError("Use pearson or spearman.")

    # Transpose so that columns are genes.
    corr = df.T.corr(method=method)

    if use_abs:
        corr = corr.abs()

    # Replace undefined correlations with zero.
    corr = corr.fillna(0.0)

    return corr


def adjacency_from_correlation(
    corr: pd.DataFrame,
    threshold: float,
    weighted: bool = False
) -> pd.DataFrame:

    if not 0 <= threshold <= 1:
        raise ValueError("Threshold must be between 0 and 1.")

    if weighted:
        adjacency = corr.where(corr >= threshold, 0.0)
    else:
        adjacency = (corr >= threshold).astype(int)

    # Remove self-connections.
    np.fill_diagonal(adjacency.values, 0)

    return adjacency


def graph_from_adjacency(
    A: pd.DataFrame,
    undirected: bool = True
) -> nx.Graph:

    if undirected:
        G = nx.from_pandas_adjacency(A)
    else:
        G = nx.from_pandas_adjacency(
            A,
            create_using=nx.DiGraph
        )

    # Remove genes without any connections.
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)

    return G


def detect_modules_louvain_or_greedy(
    G: nx.Graph
) -> Dict[str, int]:

    if G.number_of_nodes() == 0:
        return {}

    try:
        # Louvain method
        communities = nx.community.louvain_communities(
            G,
            seed=42
        )
        print("Community detection: Louvain")

    except AttributeError:
        # Fallback for older NetworkX versions
        communities = nx.community.greedy_modularity_communities(
            G
        )
        print("Community detection: Greedy modularity")

    mapping = {}

    for module_id, community in enumerate(communities):
        for gene in community:
            mapping[gene] = module_id

    return mapping


def save_modules_csv(
    mapping: Dict[str, int],
    out_csv: Path
) -> None:

    out_csv.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    df_modules = pd.DataFrame(
        {
            "Gene": list(mapping.keys()),
            "Module": list(mapping.values())
        }
    )

    df_modules = df_modules.sort_values(
        ["Module", "Gene"]
    )

    df_modules.to_csv(out_csv, index=False)


if __name__ == "__main__":

    # Load data
    df = read_expression_matrix(INPUT_CSV)
    print("Original matrix:", df.shape)

    # Normalize counts if needed
    if NORMALIZE_COUNTS:
        df = normalize_counts(df)

    # Log-transform and filter
    df_filtered = log_and_filter(
        df,
        VARIANCE_THRESHOLD
    )
    print("After filtering:", df_filtered.shape)

    # Compute correlation matrix
    corr = correlation_matrix(
        df_filtered,
        method=CORR_METHOD,
        use_abs=USE_ABS_CORR
    )

    # Create adjacency matrix
    A = adjacency_from_correlation(
        corr,
        threshold=ADJ_THRESHOLD,
        weighted=False
    )

    # Build graph
    G = graph_from_adjacency(
        A,
        undirected=MAKE_UNDIRECTED
    )

    print(
        f"Graph created with {G.number_of_nodes()} "
        f"nodes and {G.number_of_edges()} edges."
    )

    # Detect modules
    gene_to_module = detect_modules_louvain_or_greedy(G)

    print(
        f"Detected {len(set(gene_to_module.values()))} modules."
    )

    # Save output
    save_modules_csv(gene_to_module, OUTPUT_CSV)
    print(f"Saved gene-to-module mapping to: {OUTPUT_CSV}")
