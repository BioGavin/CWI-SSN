# ssn_blast.py
import subprocess
from pathlib import Path

BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send evalue bitscore qlen slen"
)

def make_db(fasta: str, db_dir: str = "results/blast") -> str:
    fasta = Path(fasta); db_dir = Path(db_dir); db_dir.mkdir(parents=True, exist_ok=True)
    db = db_dir / "seqdb"
    cmd = ["makeblastdb", "-dbtype", "prot", "-input_type", "fasta",
           "-in", str(fasta), "-parse_seqids", "-out", str(db)]
    subprocess.run(cmd, check=True, capture_output=True)
    
    return str(db)

def all_vs_all_blastp(
    fasta: str,
    db: str,
    out_file: str = "results/blast/all2all.m8",
    evalue: float = 1e5,
    threads: int = 16,
    matrix: str = "BLOSUM62",
    seg: str = "yes",
    max_hsps: int = 1,
) -> str:
    out = Path(out_file); out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["blastp", "-query", str(fasta), "-db", str(db), "-out", str(out),
           "-outfmt", BLAST_OUTFMT, "-evalue", str(evalue), "-matrix", matrix,
           "-seg", seg, "-max_hsps", str(max_hsps), "-num_threads", str(threads)]
    subprocess.run(cmd, check=True, capture_output=True)
    return str(out)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True, help="FASTA file")
    ap.add_argument("-o","--out", default="results/blast/all2all.m8")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    db = make_db(args.input)
    all_vs_all_blastp(args.input, db, args.out, threads=args.threads)
    print(f"[OK] Wrote {args.out}")