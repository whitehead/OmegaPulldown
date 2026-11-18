#!/usr/bin/env python3
"""
Generate the pair list that drives the rest of the pipeline.

A single input shape to accommodate different screen designs:

    generate_pairs.py --pairs-tsv pairs.tsv --sequences pool.fasta \
                      --output pairs.tsv

  --pairs-tsv   the pair set: an edge list whose first two columns are
                protein IDs
  --sequences   the protein pool: a FASTA holding every protein that list
                names

INPUT SCHEMA
---------------------------------------
First two columns are protein IDs, positional; a header row is detected
rather than required; further columns (MS score, FDR, spectral counts) are
ignored. Allows a raw proteomics export drop in unmodified: SAINT output
with Bait/Prey/SaintScore, or whatever a given IP-MS pipeline calls its
columns.

OUTPUT SCHEMA
----------------------
  pair_id, bait_id, bait_seq, candidate_id, candidate_seq

Pair IDs use the format {bait}_with_{candidate}  (e.g. actin_with_HuR).

CHAIN-ORDER CONTRACT
--------------------
Column 1 is the bait is chain A. The input's first column becomes bait_id
becomes chain A in the predicted complex, through every downstream stage.
The chains are physically distinguishable in the output, so something has to
pick A deterministically: this is that rule.

An AP-MS user who puts their tagged bait in column 1 therefore gets it as
chain A, which is almost always the intent.

UNORDERED PAIRS
---------------
A complex is unordered: A_with_B and B_with_A are the same system.
Pairs are collapsed on the canonical key tuple(sorted((a, b))) and the
FIRST orientation listed wins, so the collapse never overrides the
chain-order contract for the row that survives, it only discards the later
restatement.
"""

import argparse
import csv
import sys


def parse_fasta(fasta_path):
    """Yield (id, sequence) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            elif line:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def load_sequences_from_fasta(fasta_path):
    """Load all sequences from a single multi-sequence FASTA file."""
    seqs = {}
    for prot_id, seq in parse_fasta(fasta_path):
        if prot_id in seqs:
            print(f"WARNING: duplicate ID '{prot_id}' in {fasta_path}, skipping",
                  file=sys.stderr)
            continue
        seqs[prot_id] = seq
    return seqs


def read_pair_list(pairs_path, seqs):
    """
    Read (id_a, id_b) rows from a TSV or CSV. Only the first two columns are
    used; any further columns (MS score, FDR, spectral counts) are ignored
    here: they are not threaded through to the rankings yet.

    A header row is detected rather than assumed: if neither of the first
    row's two IDs resolves against the FASTA, it is treated as a header.
    That way both a bare two-column list and one exported with column names
    work without a flag.
    """
    rows = []
    with open(pairs_path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sample.splitlines()[0] else ","
        for lineno, row in enumerate(csv.reader(fh, delimiter=delim), start=1):
            row = [c.strip() for c in row if c.strip() != ""]
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                sys.exit(f"ERROR: {pairs_path} line {lineno}: need two protein "
                         f"IDs, got {row!r}")
            rows.append((row[0], row[1], lineno))

    if not rows:
        sys.exit(f"ERROR: no pairs found in {pairs_path}")

    a, b, _ = rows[0]
    if a not in seqs and b not in seqs:
        print(f"Treating first line of {pairs_path} as a header: {a!r}, {b!r}",
              file=sys.stderr)
        rows = rows[1:]
        if not rows:
            sys.exit(f"ERROR: {pairs_path} contains only a header")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pairs-tsv", dest="pairs_tsv", required=True,
                        help="TSV/CSV whose first two columns name the pairs "
                             "to predict. Column 1 becomes chain A.")
    parser.add_argument("--sequences", required=True,
                        help="FASTA holding every protein referenced by --pairs-tsv")
    parser.add_argument("--output", default="pairs.tsv", help="Output TSV file")
    args = parser.parse_args()

    seqs = load_sequences_from_fasta(args.sequences)
    if not seqs:
        sys.exit(f"ERROR: no FASTA sequences found in {args.sequences}")
    rows = read_pair_list(args.pairs_tsv, seqs)
    print(f"Loaded {len(seqs)} sequences and {len(rows)} listed pairs",
          file=sys.stderr)

    # Resolve every ID before writing anything. A pair list naming a protein
    # the pool does not hold is a specification error.
    missing = {}
    for a, b, lineno in rows:
        for pid in (a, b):
            if pid not in seqs:
                missing.setdefault(pid, lineno)
    if missing:
        shown = list(missing.items())[:10]
        detail = "\n".join(f"    {pid}  (first seen line {ln})"
                           for pid, ln in shown)
        more = "" if len(missing) <= 10 else f"\n    ... and {len(missing) - 10} more"
        sys.exit(f"ERROR: {len(missing)} protein ID(s) in {args.pairs_tsv} "
                 f"are absent from {args.sequences}:\n{detail}{more}")

    seen = set()
    dupes = 0
    n = 0
    with open(args.output, "w") as out:
        out.write("pair_id\tbait_id\tbait_seq\tcandidate_id\tcandidate_seq\n")
        for a, b, _ in rows:
            key = tuple(sorted((a, b)))
            if key in seen:
                dupes += 1
                continue
            seen.add(key)
            out.write(f"{a}_with_{b}\t{a}\t{seqs[a]}\t{b}\t{seqs[b]}\n")
            n += 1

    # Proteins in the pool that no surviving pair names get no MSA, because
    # the MSA universe is read from this file (main.nf) rather than from the
    # pool.
    used = {pid for pair in seen for pid in pair}
    unused = len(seqs) - len(used)
    if unused:
        print(f"{unused} of {len(seqs)} pooled protein(s) are named by no pair; "
              f"they get no MSA and no prediction", file=sys.stderr)
    if dupes:
        print(f"Collapsed {dupes} duplicate unordered pair(s) "
              f"(A_with_B == B_with_A)", file=sys.stderr)
    print(f"Wrote {n} pairs to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
