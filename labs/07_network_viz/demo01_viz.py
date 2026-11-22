"""
Demo 1 — Visualization of a toy network with modules and hub genes

Goal:
- Show how to color nodes based on module
- How to identify hub genes (high degree)
- How to export the figure

This demo uses a small toy graph, not real data.
"""

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# 1) Build a small graph
G = nx.Graph()
edges = [
    ("GeneA", "GeneB"),
    ("GeneA", "GeneC"),
    ("GeneB", "GeneC"),
    ("GeneC", "GeneD"),
    ("GeneD", "GeneE"),
    ("GeneE", "GeneF"),
]
G.add_edges_from(edges)

# 2) Modules assigned manually
gene2module = {
    "GeneA": 1, "GeneB": 1, "GeneC": 1,
    "GeneD": 2, "GeneE": 2, "GeneF": 2,
}

# 3) Colors by module
cmap = plt.get_cmap("tab10")
node_colors = [cmap((gene2module[n] - 1) % 10) for n in G.nodes()]

# 4) Hub genes (high degree)
deg = dict(G.degree())
hubs = sorted(deg.items(), key=lambda x: x[1], reverse=True)[:2]  # top 2
hub_nodes = {n for n, _ in hubs}
print("Hub genes:", hubs)

# 5) Layout + visualization
pos = nx.spring_layout(G, seed=42)
plt.figure(figsize=(6, 5))
nx.draw_networkx_edges(G, pos, alpha=0.3)
nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    node_size=[300 if n in hub_nodes else 150 for n in G.nodes()]
)
nx.draw_networkx_labels(G, pos, font_size=10)

plt.axis("off")
plt.tight_layout()
plt.savefig("labs/07_network_viz/demo_network.png", dpi=200)
plt.show()

# 6) Export hub genes to CSV
df_hubs = pd.DataFrame(hubs, columns=["Gene", "Degree"])
df_hubs.to_csv("labs/07_network_viz/demo_hubs.csv", index=False)
print("Saved demo_network.png and demo_hubs.csv")
