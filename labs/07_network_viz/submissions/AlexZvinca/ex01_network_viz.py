
"""
Week 7 — Visualization of Gene Co-Expression Networks

Reconstruct the Week 6 network, color genes by module,
highlight hub genes, and save the network visualization.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import networkx as nx

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------
# Configuration
# --------------------------
HANDLE = "AlexZvinca"

EXPR_CSV = Path(
    f"data/work/{HANDLE}/lab06/expression_matrix.csv"
)

MODULES_CSV = Path(
    f"labs/07_network_viz/modules_{HANDLE}.csv"
)

PRECOMPUTED_ADJ_CSV: Optional[Path] = None

# Same settings used in Week 6
CORR_METHOD = "spearman"
USE_ABS_CORR = True
VARIANCE_THRESHOLD = 0.5
ADJ_THRESHOLD = 0.7
WEIGHTED = False

# Visualization
SEED = 42
TOPK_HUBS = 10
NODE_BASE_SIZE = 300
EDGE_ALPHA = 0.25

OUT_DIR = Path(
    f"labs/07_network_viz/submissions/{HANDLE}"
)

OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_HUBS = OUT_DIR / f"hubs_{HANDLE}.csv"


# --------------------------
# 1. Check input files
# --------------------------
def ensure_exists(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")


# --------------------------
# 2. Read expression matrix
# --------------------------
def read_expression_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)

    df = df.apply(pd.to_numeric, errors="raise")

    if df.empty:
        raise ValueError("Expression matrix is empty.")

    if df.index.has_duplicates:
        raise ValueError("Duplicate gene IDs found.")

    if not np.isfinite(df.to_numpy()).all():
        raise ValueError("Invalid expression values found.")

    if (df < 0).any().any():
        raise ValueError("Negative expression values found.")

    return df


# --------------------------
# 3. Read module assignments
# --------------------------
def read_modules_csv(path: Path) -> Dict[str, int]:
    df = pd.read_csv(path)

    if not {"Gene", "Module"}.issubset(df.columns):
        raise ValueError("CSV must contain Gene and Module columns.")

    if df["Gene"].duplicated().any():
        raise ValueError("Duplicate genes in module CSV.")

    return dict(zip(df["Gene"], df["Module"]))


# --------------------------
# 4. Preprocess expression
# --------------------------
def preprocess_expression(df: pd.DataFrame) -> pd.DataFrame:
    # Our Week 6 input was already processed expression data.
    # We therefore do not apply CPM normalization again.

    df_log = np.log2(df + 1)

    variances = df_log.var(axis=1)

    df_filtered = df_log.loc[
        variances > VARIANCE_THRESHOLD
    ]

    if len(df_filtered) < 2:
        raise ValueError("Too few genes remain after filtering.")

    return df_filtered


# --------------------------
# 5. Correlation to adjacency
# --------------------------
def correlation_to_adjacency(
    expr: pd.DataFrame,
    method: str,
    use_abs: bool,
    threshold: float,
    weighted: bool
) -> pd.DataFrame:

    corr = expr.T.corr(method=method)

    if use_abs:
        corr = corr.abs()

    corr = corr.fillna(0.0)

    if weighted:
        A = corr.where(corr >= threshold, 0.0)
    else:
        A = (corr >= threshold).astype(int)

    np.fill_diagonal(A.values, 0)

    return A


# --------------------------
# 6. Build network graph
# --------------------------
def graph_from_adjacency(A: pd.DataFrame) -> nx.Graph:
    G = nx.from_pandas_adjacency(A)

    # Week 6 removed isolated genes.
    G.remove_nodes_from(list(nx.isolates(G)))

    return G


# --------------------------
# 7. Module colors
# --------------------------
def color_map_from_modules(
    nodes: Iterable[str],
    gene2module: Dict[str, int]
) -> Dict[str, tuple]:

    cmap = plt.get_cmap("tab20")

    return {
        gene: cmap(int(gene2module[gene]) % 20)
        for gene in nodes
    }


# --------------------------
# 8. Identify hub genes
# --------------------------
def compute_hubs(
    G: nx.Graph,
    topk: int
) -> pd.DataFrame:

    degrees = dict(G.degree())

    hubs = sorted(
        degrees.items(),
        key=lambda item: (-item[1], item[0])
    )[:topk]

    return pd.DataFrame(
        hubs,
        columns=["Gene", "Degree"]
    )


# --------------------------
# Main
# --------------------------
def main():
    # Verify inputs
    ensure_exists(EXPR_CSV)
    ensure_exists(MODULES_CSV)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load Week 6 data
    expr = read_expression_matrix(EXPR_CSV)
    gene2module = read_modules_csv(MODULES_CSV)

    print("Original expression matrix:", expr.shape)
    print("Genes with module assignments:", len(gene2module))

    # Reproduce Week 6 preprocessing
    expr = preprocess_expression(expr)
    print("After variance filtering:", expr.shape)

    # Load or reconstruct adjacency
    if PRECOMPUTED_ADJ_CSV is not None:
        ensure_exists(PRECOMPUTED_ADJ_CSV)

        A = pd.read_csv(
            PRECOMPUTED_ADJ_CSV,
            index_col=0
        )
    else:
        A = correlation_to_adjacency(
            expr,
            method=CORR_METHOD,
            use_abs=USE_ABS_CORR,
            threshold=ADJ_THRESHOLD,
            weighted=WEIGHTED
        )

    # Construct the full graph first
    G = graph_from_adjacency(A)

    # Verify the reconstructed graph matches Week 6 modules
    graph_genes = set(G.nodes())
    module_genes = set(gene2module)

    if graph_genes != module_genes:
        missing = module_genes - graph_genes
        extra = graph_genes - module_genes

        raise ValueError(
            "Week 7 graph does not match Week 6 modules.\n"
            f"Genes missing from graph: {sorted(missing)}\n"
            f"Extra genes in graph: {sorted(extra)}\n"
            "Check that preprocessing and thresholds match Week 6."
        )

    print(
        f"Graph created with {G.number_of_nodes()} nodes "
        f"and {G.number_of_edges()} edges."
    )

    print("Modules:", len(set(gene2module.values())))

    # Colors
    color_map = color_map_from_modules(
        G.nodes(),
        gene2module
    )

    node_colors = [
        color_map[gene]
        for gene in G.nodes()
    ]

    # Hub genes
    hubs_df = compute_hubs(G, TOPK_HUBS)
    hub_nodes = set(hubs_df["Gene"])

    print("\nTop hub genes:")
    print(hubs_df.to_string(index=False))

    node_sizes = [
        NODE_BASE_SIZE * 2.5
        if gene in hub_nodes
        else NODE_BASE_SIZE
        for gene in G.nodes()
    ]

    # Spring layout
    pos = nx.spring_layout(G, seed=SEED)

    # Draw graph
    fig, ax = plt.subplots(figsize=(13, 10))

    nx.draw_networkx_edges(
        G,
        pos,
        alpha=EDGE_ALPHA,
        width=1.3,
        ax=ax
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        edgecolors="black",
        linewidths=0.5,
        ax=ax
    )

    # Label only hub genes
    hub_labels = {
        gene: gene
        for gene in hub_nodes
    }

    nx.draw_networkx_labels(
        G,
        pos,
        labels=hub_labels,
        font_size=9,
        font_weight="bold",
        ax=ax
    )

    ax.set_title(
        "Gene Co-Expression Network — Week 7\n"
        "Node color = module | Larger nodes = hub genes",
        fontsize=14
    )

    ax.axis("off")

    plt.tight_layout()
    plt.savefig(
        OUT_PNG,
        dpi=200,
        bbox_inches="tight"
    )
    plt.close(fig)

    # Save hub genes
    hubs_df.to_csv(OUT_HUBS, index=False)

    print("\nSaved network figure:", OUT_PNG)
    print("Saved hub gene table:", OUT_HUBS)


if __name__ == "__main__":
    main()
