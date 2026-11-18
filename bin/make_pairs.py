#!/usr/bin/env python3
"""
Evaluate a pair-generating RULE into an explicit pair list.

This is a PRE-PROCESSING tool, not a pipeline stage. It runs before
`nextflow run`, and its output is fed back in as --pairs_tsv.

    make_pairs.py --all-vs-all pool.fasta            > pairs.tsv
    make_pairs.py --baits one.txt --candidates c.txt > pairs.tsv

OUTPUT
------
A two-column TSV -- the INPUT pair-list schema, i.e. an edge list:

    bait_id<TAB>candidate_id

Note this is NOT the same file as the pipeline's emitted
`results/pairs/pairs.tsv`, which additionally inlines both sequences. This
file is an input; that one is the initialized screen. GENERATE_PAIRS is the
resolver between them.

CHAIN-ORDER
--------------------
Column 1 is the bait and defined as chain A. That holds through
generate_pairs.py and into every prediction, so the orientation this tool
writes is the orientation delivered in the output structures.

Both rules below are deterministic:

  --all-vs-all           column 1 is the lexicographically smaller ID
  --baits/--candidates   column 1 is the protein from the baits file

Pairs are unordered complexes, so A_with_B and B_with_A are collapsed on the
canonical key tuple(sorted((a, b))) and the first orientation emitted wins.
With disjoint bait and candidate sets that is a no-op; with overlapping sets
it is exactly the N^2 -> N(N-1)/2 halving.
"""

import argparse
import sys
from itertools import combinations, product


def parse_fasta_ids(path):
    """Yield protein IDs from a FASTA file, in file order.

    The ID is the first whitespace-delimited token of the header, matching
    bin/generate_pairs.py's parse_fasta(): the two must agree or a protein
    named here will not resolve against the pool later.
    """
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def load_ids(path):
    """Load protein IDs from either a FASTA or a plain one-ID-per-line list.

    Determined by content rather than by file extension, so a bait set can be
    given as a FASTA on hand or as a bare list of IDs from a spreadsheet,
    whichever is closer to hand. Blank lines and #-comments are skipped in
    list form.
    """
    with open(path) as fh:
        lines = fh.readlines()

    first = next((l for l in lines if l.strip()), None)
    if first is None:
        sys.exit(f"ERROR: {path} is empty")

    if first.startswith(">"):
        return [l[1:].split()[0] for l in lines if l.startswith(">")]

    ids = []
    for line in lines:
        # Take the first field, so a two-column list (ID + anything) also
        # works without pre-trimming.
        tok = line.split("#", 1)[0].split()
        if tok:
            ids.append(tok[0])
    return ids


def dedup_preserving_order(ids, label):
    """Drop repeated IDs, keeping first occurrence; warn once per duplicate."""
    seen = set()
    out = []
    for i in ids:
        if i in seen:
            print(f"WARNING: duplicate ID '{i}' in {label}, ignoring the repeat",
                  file=sys.stderr)
            continue
        seen.add(i)
        out.append(i)
    return out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--all-vs-all", dest="all_vs_all", metavar="POOL.fasta",
                        help="Pair every protein in this FASTA with every other "
                             "-- N(N-1)/2 pairs.")
    parser.add_argument("--baits", metavar="FILE",
                        help="Bait set, as a FASTA or a one-ID-per-line list "
                             "(sniffed). Each bait is paired with every candidate.")
    parser.add_argument("--candidates", metavar="FILE",
                        help="Candidate set, as a FASTA or a one-ID-per-line list.")
    parser.add_argument("--sequences", metavar="POOL.fasta",
                        help="Optional with --baits/--candidates: the pool the "
                             "pipeline will be run against. Every ID is checked "
                             "against it now, so a typo fails here in a second "
                             "rather than at GENERATE_PAIRS.")
    parser.add_argument("--self-pairs", action="store_true",
                        help="Also emit homodimers (A with A). Off by default.")
    parser.add_argument("--output", default="-",
                        help="Output TSV ('-' for stdout, the default).")
    args = parser.parse_args()

    # --- Rule selection: exactly one rule, stated exactly once -------------
    if args.all_vs_all and (args.baits or args.candidates):
        sys.exit("ERROR: --all-vs-all is mutually exclusive with "
                 "--baits/--candidates (they are two different rules).")
    if not args.all_vs_all and not (args.baits and args.candidates):
        sys.exit("ERROR: need either --all-vs-all POOL.fasta, or both --baits "
                 "and --candidates.")

    pool = None
    if args.sequences:
        pool = set(parse_fasta_ids(args.sequences))
        if not pool:
            sys.exit(f"ERROR: no FASTA sequences found in {args.sequences}")

    if args.all_vs_all:
        ids = dedup_preserving_order(parse_fasta_ids(args.all_vs_all),
                                     args.all_vs_all)
        if not ids:
            sys.exit(f"ERROR: no FASTA sequences found in {args.all_vs_all}")
        if len(ids) < 2 and not args.self_pairs:
            sys.exit(f"ERROR: {args.all_vs_all} holds {len(ids)} protein(s); "
                     f"an all-vs-all needs at least 2 (or --self-pairs).")
        # Lexicographic order IS the chain-order rule here: sorting first
        # means combinations() emits (smaller, larger), so the smaller ID is
        # column 1 is chain A. Deterministic, and independent of the order
        # the FASTA happens to be in.
        ids = sorted(ids)
        rule_pairs = combinations(ids, 2)
        if args.self_pairs:
            rule_pairs = list(rule_pairs) + [(i, i) for i in ids]
        described = f"all-vs-all over {len(ids)} proteins from {args.all_vs_all}"
    else:
        baits = dedup_preserving_order(load_ids(args.baits), args.baits)
        cands = dedup_preserving_order(load_ids(args.candidates), args.candidates)
        if not baits:
            sys.exit(f"ERROR: no protein IDs found in {args.baits}")
        if not cands:
            sys.exit(f"ERROR: no protein IDs found in {args.candidates}")
        # Baits iterate in the OUTER loop, so a bait reaches the canonical key
        # first and therefore owns column 1: the bait is chain A, which is
        # the point of calling it the bait. Both sides are sorted so the
        # output is reproducible regardless of input file ordering.
        rule_pairs = product(sorted(baits), sorted(cands))
        described = (f"{len(baits)} baits x {len(cands)} candidates "
                     f"({args.baits} x {args.candidates})")

    # --- Validate against the pool, if one was given -----------------------
    if pool is not None:
        named = set()
        if args.all_vs_all:
            named.update(ids)
        else:
            named.update(baits)
            named.update(cands)
        missing = sorted(named - pool)
        if missing:
            shown = "\n".join(f"    {m}" for m in missing[:10])
            more = "" if len(missing) <= 10 else f"\n    ... and {len(missing) - 10} more"
            sys.exit(f"ERROR: {len(missing)} protein ID(s) are absent from "
                     f"{args.sequences}:\n{shown}{more}")

    # --- Emit --------------------------------------------------------------
    out = sys.stdout if args.output == "-" else open(args.output, "w")
    seen = set()
    dupes = 0
    n = 0
    try:
        out.write("bait_id\tcandidate_id\n")
        for a, b in rule_pairs:
            if a == b and not args.self_pairs:
                continue
            key = tuple(sorted((a, b)))
            if key in seen:
                dupes += 1
                continue
            seen.add(key)
            out.write(f"{a}\t{b}\n")
            n += 1
    finally:
        if out is not sys.stdout:
            out.close()

    if dupes:
        print(f"Collapsed {dupes} duplicate unordered pair(s) "
              f"(A_with_B == B_with_A)", file=sys.stderr)
    dest = "stdout" if args.output == "-" else args.output
    print(f"Wrote {n} pairs to {dest}  [{described}]", file=sys.stderr)


if __name__ == "__main__":
    main()
