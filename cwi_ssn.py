#!/usr/bin/env python3
"""Command-line driver that stitches the modified SSN pipeline together."""

import argparse
from pathlib import Path

from ssn_blastp import make_db, all_vs_all_blastp
from ssn_build import reduce_edges, add_scores, build_graph


def load_sequences(fasta_path: Path) -> dict:
    seqs = {}
    current_id = None
    fragments = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(fragments)
                current_id = line[1:].split()[0]
                fragments = []
            else:
                if current_id is None:
                    raise ValueError("FASTA format error: sequence line before header")
                fragments.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(fragments)
    return seqs


def run_cwi_ssn(
    fasta: str,
    output_root: str = "results",
    threads: int = 16,
    cwi_threshold: float = 0.6,
    include_singletons: bool = False,
) -> dict:
    """Run the entire CWI-SSN construction pipeline and return generated paths."""
    fasta = Path(fasta)
    workdir = Path(output_root)
    blast_dir = workdir / "blast"
    ssn_dir = workdir / "ssn"
    blast_dir.mkdir(parents=True, exist_ok=True)
    ssn_dir.mkdir(parents=True, exist_ok=True)

    threshold_tag = f"{cwi_threshold}".rstrip("0").rstrip(".") or "0"
    graphml_path = ssn_dir / f"network_cwi{threshold_tag}.graphml"
    evalue = 1e5
    matrix = "BLOSUM62"
    sequences = load_sequences(fasta)
    db_path = make_db(str(fasta), db_dir=str(blast_dir))
    print(f"[1/5] BLAST database ready: {db_path}")

    m8_path = all_vs_all_blastp(
        str(fasta),
        db_path,
        out_file=str(blast_dir / "all2all.m8"),
        evalue=evalue,
        matrix=matrix,
        threads=threads,
    )
    print(f"[2/5] All-vs-all BLAST finished: {m8_path}")

    reduced_path = reduce_edges(m8_path, str(blast_dir / "reduced.tsv"))
    print(f"[3/5] Reduced to one edge per sequence pair: {reduced_path}")

    scored_path = add_scores(reduced_path, str(ssn_dir / "edges_scored.tsv"))
    print(f"[4/5] Added coverage-weighted identity scores: {scored_path}")

    graph_path = build_graph(
        scored_path,
        cwi_threshold,
        str(graphml_path),
        node_sequences=sequences,
        include_singletons=include_singletons,
    )
    print(f"[5/5] GraphML network written: {graph_path}")

    return {
        "db": db_path,
        "blast_m8": m8_path,
        "reduced_edges": reduced_path,
        "scored_edges": scored_path,
        "graphml": graph_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Construct a coverage-weighted identity (CWI) sequence similarity "
            "network starting from a protein FASTA file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Protein FASTA used both to build the BLAST database and as BLAST queries.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        "--output-root",
        dest="output_root",
        default="results",
        help="Base directory where intermediate BLAST outputs and SSN files are stored.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of CPU threads passed to BLASTP.",
    )
    parser.add_argument(
        "--cwi-threshold",
        type=float,
        default=0.6,
        help="Minimum coverage-weighted identity for keeping an edge in the SSN.",
    )
    parser.add_argument(
        "--include-singletons",
        action="store_true",
        help="When set, retain sequences without edges as singleton nodes.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_cwi_ssn(
        fasta=args.input,
        output_root=args.output_root,
        threads=args.threads,
        cwi_threshold=args.cwi_threshold,
        include_singletons=args.include_singletons,
    )
    print("[OK] Finished building CWI SSN:")
    for label, path in outputs.items():
        print(f"  {label:>14}: {path}")


if __name__ == "__main__":
    main()
