#!/usr/bin/env python3
import argparse, re
from io import StringIO
from typing import Set, List
import pandas as pd
from Bio import Phylo

def read_tree(path):
    txt = open(path, "r").read()
    # Strip FigTree/BEAST attributes
    txt = re.sub(r"\[\&[^\]]*\]", "", txt)
    fmt = "nexus" if "#NEXUS" in txt[:200].upper() else "newick"
    return Phylo.read(StringIO(txt), fmt)

def tip_names(clade) -> Set[str]:
    return set(t.name for t in clade.get_terminals() if t.name)

def root_split_from_timetree(ttree):
    """Return the two tip-sets below the time-tree root’s two children."""
    root = ttree.root
    if len(root.clades) < 2:
        raise ValueError("Time tree root has <2 children; cannot define a split.")
    # Choose the two largest child clades if more than 2
    children = sorted(root.clades, key=lambda c: len(list(c.get_terminals())), reverse=True)
    A = tip_names(children[0])
    B = tip_names(children[1])
    return A, B

def mrca_of_names(tree, names: Set[str]):
    # Take only names present in this tree
    present = [t for t in tree.get_terminals() if t.name in names]
    if not present:
        raise ValueError("None of the provided tips are present to define the MRCA.")
    return tree.common_ancestor(present)

def reroot_divergence_to_match(div_tree, ttime_tree):
    A, B = root_split_from_timetree(ttime_tree)
    # Pick the smaller side as an "outgroup" to induce the same split
    outgroup_set = A if len(A) <= len(B) else B
    # If nothing overlaps (name differences), try the other side
    try:
        og_mrca = mrca_of_names(div_tree, outgroup_set)
    except Exception:
        other = B if outgroup_set is A else A
        og_mrca = mrca_of_names(div_tree, other)
    # Root with that clade as outgroup (places root between the two split sides)
    div_tree.root_with_outgroup(og_mrca)
    return div_tree

def root_to_tip_df(tree):
    if not tree.rooted:
        raise ValueError("Tree is not rooted.")
    root = tree.root
    data = []
    for tip in tree.get_terminals():
        d = tree.distance(root, tip)
        d = float(d) if d is not None else 0.0
        data.append((tip.name, d))
    return pd.DataFrame(data, columns=["name", "root_to_tip"])

def detect_cols(df, name_col=None, date_col=None):
    df.columns = [c.strip() for c in df.columns]
    if name_col is None:
        for c in ["name","strain","seqName","taxon","Virus name"]:
            if c in df.columns: name_col=c; break
    if date_col is None:
        for c in ["date_decimal","decimal_date","decimal-year","decimal-date","date"]:
            if c in df.columns: date_col=c; break
    if name_col is None or date_col is None:
        raise ValueError(f"Could not infer columns from {list(df.columns)}; "
                         f"use --name-column / --date-column.")
    return name_col, date_col

def main():
    ap = argparse.ArgumentParser(description="Recreate RTT table exactly from TreeTime outputs.")
    ap.add_argument("--divergence-tree", required=True,
                    help="treetime_output/divergence_tree.nexus")
    ap.add_argument("--timetree", required=True,
                    help="treetime_output/timetree.nexus (defines the root split)")
    ap.add_argument("--dates", required=True,
                    help="metadata/dates.tsv (decimal years or ISO)")
    ap.add_argument("--out", default="rtt_from_fullrun.csv")
    ap.add_argument("--name-column", default=None)
    ap.add_argument("--date-column", default=None)
    ap.add_argument("--skiprows", type=int, default=0)
    args = ap.parse_args()

    # 1) read trees
    div_tree = read_tree(args.divergence_tree)
    tt_tree  = read_tree(args.timetree)

    # 2) re-root divergence tree to match time-tree root split
    div_tree = reroot_divergence_to_match(div_tree, tt_tree)

    # 3) compute root-to-tip (genetic)
    rtt = root_to_tip_df(div_tree)

    # 4) merge dates
    dates = pd.read_csv(args.dates, sep="\t", dtype=str, skiprows=args.skiprows).fillna("")
    name_col, date_col = detect_cols(dates, args.name_column, args.date_column)
    # try numeric date
    dates[date_col] = pd.to_numeric(dates[date_col], errors="coerce")
    rtt["name"] = rtt["name"].astype(str).str.strip()
    dates[name_col] = dates[name_col].astype(str).str.strip()
    out = rtt.merge(dates[[name_col, date_col]], left_on="name", right_on=name_col, how="left")
    out = out.rename(columns={date_col: "date"}).drop(columns=[name_col])

    out.to_csv(args.out, index=False)
    print(f"Wrote {len(out)} rows → {args.out}")
    missing = out["date"].isna().sum()
    if missing:
        print(f"Note: {missing} tips have no date in the TSV (left as NaN).")

if __name__ == "__main__":
    main()