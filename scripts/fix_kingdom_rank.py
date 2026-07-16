#!/usr/bin/env python3
"""
Backfill the missing kingdom/domain-level entry in a Taxor taxonomy input CSV.

Since NCBI stopped ranking Bacteria/Archaea/Eukaryota as "superkingdom" (now
"domain") and Viruses as something else entirely ("acellular root"), tools
like `taxonkit reformat` that build the k__/p__/c__/.../s__ lineage string by
matching the *rank name* "superkingdom" emit an empty (or entirely missing)
kingdom slot for every organism.

This re-derives the kingdom/domain level directly from taxdump's parent-child
structure (nodes.dmp), which is unaffected by rank renaming, using a small
set of well-known, stable anchor taxids, and patches the taxnames_string /
taxid_string columns of the Taxor input CSV so they always start with a
correctly populated "k__<Name>" entry.

Usage:
    python3 fix_kingdom_rank.py \
        --nodes taxdump/nodes.dmp \
        --input refseq_accessions_taxonomy.csv \
        --output refseq_accessions_taxonomy.fixed.csv

Assumes the Taxor input CSV is tab-separated with columns (0-indexed):
  0 accession_id, 1 taxid, 2 path, 3 organism_name,
  4 taxnames_string (k__..;p__..;...;s__..), 5 taxid_string (taxid;taxid;...)
matching the format documented in Taxor's README.
"""
import argparse
import os
import sys

# Stable, rank-name-independent anchor taxids for the domain/kingdom level.
ANCHOR_TAXIDS = {
    "2": "Bacteria",
    "2157": "Archaea",
    "2759": "Eukaryota",
    "10239": "Viruses",
    "12884": "Viroids",
    "12908": "unclassified sequences",
    "28384": "other sequences",
}

RANK_PREFIXES = ["k", "p", "c", "o", "f", "g", "s"]  # canonical 7-level order


def load_parents(nodes_path):
    parent = {}
    with open(nodes_path) as f:
        for line in f:
            fields = line.split("\t|\t")
            taxid = fields[0].strip()
            parent_taxid = fields[1].strip()
            parent[taxid] = parent_taxid
    return parent


def find_domain_anchor(taxid, parent):
    """Walk up the parent chain from taxid to root; return the first anchor
    taxid encountered, or None if none of the known anchors are ancestors."""
    seen = set()
    current = taxid
    while current and current not in seen:
        if current in ANCHOR_TAXIDS:
            return current
        seen.add(current)
        nxt = parent.get(current)
        if nxt is None or nxt == current:  # reached root (root's parent is itself)
            break
        current = nxt
    return None


def split_ranked_string(value):
    """Split a 'k__X;p__Y;...' string into {prefix: name}, keyed by each
    segment's own k__/p__/.../s__ marker (not position) - a rank slot may be
    fully dropped rather than left present-but-empty, and we don't want to
    assume which."""
    result = {}
    if not value:
        return result
    for segment in value.split(";"):
        if len(segment) >= 2 and segment[1] == "_":
            result[segment[0]] = segment[3:] if len(segment) > 3 else ""
    return result


def fix_row(taxid, name_string, taxid_string, parent):
    anchor = find_domain_anchor(taxid, parent)
    if anchor is None:
        # No known domain anchor found in this taxid's ancestry - leave as-is
        # rather than guess.
        return name_string, taxid_string, False

    names_by_rank = split_ranked_string(name_string)
    # taxid_string has no per-segment markers, so pair it up positionally
    # against whichever ranks the name string says are present, in canonical
    # order.
    present_ranks_in_order = [p for p in RANK_PREFIXES if p in names_by_rank]
    taxid_values = taxid_string.split(";") if taxid_string else []
    taxids_by_rank = dict(zip(present_ranks_in_order, taxid_values))

    changed = names_by_rank.get("k", "") == "" or taxids_by_rank.get("k", "") == ""

    names_by_rank["k"] = ANCHOR_TAXIDS[anchor]
    taxids_by_rank["k"] = anchor

    new_name_string = ";".join(f"{p}__{names_by_rank.get(p, '')}" for p in RANK_PREFIXES)
    new_taxid_string = ";".join(str(taxids_by_rank.get(p, "")) for p in RANK_PREFIXES)
    return new_name_string, new_taxid_string, changed


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nodes", required=True, help="path to taxdump/nodes.dmp")
    ap.add_argument("--input", required=True, help="existing Taxor taxonomy input CSV")
    ap.add_argument("--output", required=True, help="path to write the patched CSV")
    args = ap.parse_args()

    if os.path.abspath(args.input) == os.path.abspath(args.output):
        # opening --output for writing truncates it immediately, which would wipe
        # --input out from under the still-open read handle if they're the same file
        sys.exit("error: --input and --output must be different files "
                 "(write to a temp file and mv it into place afterwards)")

    parent = load_parents(args.nodes)

    fixed = 0
    total = 0
    no_anchor = 0
    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            total += 1
            cols = line.split("\t")
            if len(cols) < 6:
                fout.write(line + "\n")
                continue
            taxid = cols[1]
            name_string, taxid_string, changed = fix_row(taxid, cols[4], cols[5], parent)
            if changed:
                fixed += 1
            if name_string == cols[4] and taxid_string == cols[5] and find_domain_anchor(taxid, parent) is None:
                no_anchor += 1
            cols[4] = name_string
            cols[5] = taxid_string
            fout.write("\t".join(cols) + "\n")

    print(f"Processed {total} rows: patched kingdom-level entry for {fixed}; "
          f"{no_anchor} had no recognized domain anchor in their ancestry (left unchanged).",
          file=sys.stderr)


if __name__ == "__main__":
    main()
