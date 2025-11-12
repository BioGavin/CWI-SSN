# ssn_blast.py
from pathlib import Path
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastpCommandline

BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send evalue bitscore qlen slen"
)

def make_db(fasta: str, db_dir: str = "results/blast") -> str:
    fasta = Path(fasta); db_dir = Path(db_dir); db_dir.mkdir(parents=True, exist_ok=True)
    db = db_dir / "seqdb"
    NcbimakeblastdbCommandline(dbtype="prot", input_file=str(fasta),
                               parse_seqids=True, out=str(db))()
    
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
    cl = NcbiblastpCommandline(
        query=str(fasta), db=str(db), out=str(out),
        outfmt=BLAST_OUTFMT, evalue=evalue, matrix=matrix,
        seg=seg, max_hsps=max_hsps, num_threads=threads
    )
    cl()[0]
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