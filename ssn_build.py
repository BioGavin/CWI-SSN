# ssn_build.py
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

COLS = ["q","s","pident","alnlen","mismatch","gapopen","qstart","qend",
        "sstart","send","e","bits","qlen","slen"]

def reduce_edges(m8_path: str, out_path: str) -> str:
    df = pd.read_csv(m8_path, sep="\t", header=None, names=COLS)
    df = df[df.q != df.s].copy()
    df["u"] = df[["q","s"]].min(axis=1)
    df["v"] = df[["q","s"]].max(axis=1)
    # Pick the best hit per unordered pair using lower e-value / higher bitscore.
    red = (df.sort_values(["u","v","e","bits"], ascending=[True,True,True,False])
             .groupby(["u","v"], as_index=False).first())
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    red.to_csv(out_path, sep="\t", index=False)
    return out_path

def add_scores(edge_path: str, out_path: str) -> str:
    df = pd.read_csv(edge_path, sep="\t")
    cov_geo = np.sqrt((df["alnlen"] / df["qlen"]) * (df["alnlen"] / df["slen"]))
    ident = df["pident"] / 100.0
    df["cov_geo"] = cov_geo
    df["CWI"] = ident * cov_geo
    df.to_csv(out_path, sep="\t", index=False)
    return out_path

def build_graph(
    edge_tsv: str,
    cwi_threshold: float,
    out_graphml: str,
    node_sequences: dict | None = None,
    include_singletons: bool = True,
) -> str:
    df = pd.read_csv(edge_tsv, sep="\t")

    # Filter by the requested CWI threshold.
    df = df[df["CWI"] >= cwi_threshold].copy()

    G = nx.Graph()
    for _, r in df.iterrows():
        G.add_edge(
            r["u"], r["v"],
            CWI=float(r["CWI"]),
            cov_geo=float(r.get("cov_geo", np.nan)),
            pident=float(r["pident"]),
            bits=float(r["bits"]),
            alnlen=float(r["alnlen"]),
            qlen=float(r["qlen"]),
            slen=float(r["slen"])
        )

    if node_sequences:
        for node_id, seq in node_sequences.items():
            if include_singletons:
                G.add_node(node_id)
            if node_id in G and seq is not None:
                G.nodes[node_id]["sequence"] = seq

    Path(out_graphml).parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(G, out_graphml)
    return out_graphml

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--m8", default="results/blast/all2all.m8")
    ap.add_argument("--reduced", default="results/blast/reduced.tsv")
    ap.add_argument("--scored", default="results/ssn/edges_scored.tsv")
    ap.add_argument("--cwi_threshold", type=float, default=0.6)  # Only uses the CWI cutoff.
    ap.add_argument("--graphml", default="results/ssn/network.graphml")
    args = ap.parse_args()

    reduce_edges(args.m8, args.reduced)
    add_scores(args.reduced, args.scored)
    build_graph(args.scored, args.cwi_threshold, args.graphml)
    print(f"[OK] Wrote {args.graphml}")
