import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re

# ── 1. Read BLASTn output ────────────────────────────────────────────────────
cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
blast = pd.read_csv("./pvir_bovis.tsv", sep="\t", header=None, names=cols)

def clean_name(x):
    x = x.replace("NC_003277.2_plasmid_", "")
    x = re.sub(r"_contig_\d+", "", x)
    return x

blast["query"]   = blast["qseqid"].apply(clean_name)
blast["subject"] = blast["sseqid"].apply(clean_name)

all_plasmids = sorted(set(blast["query"].tolist() + blast["subject"].tolist()))

# ── 2. Load metadata and validate ───────────────────────────────────────────
meta = pd.read_excel("./IncFIB_plasmids_information.xlsx")
meta["label"] = meta["FILE"].apply(clean_name)

in_blast_not_meta = set(all_plasmids) - set(meta["label"])
in_meta_not_blast = set(meta["label"]) - set(all_plasmids)

print(f"Sequences in BLAST output:  {len(all_plasmids)}")
print(f"Sequences in metadata:      {len(meta)}")

if in_blast_not_meta:
    print(f"\n⚠️  In BLAST but MISSING from metadata ({len(in_blast_not_meta)}):")
    for s in sorted(in_blast_not_meta):
        print(f"   {s}")
else:
    print("\n✅  All BLAST sequences found in metadata.")

if in_meta_not_blast:
    print(f"\n⚠️  In metadata but MISSING from BLAST ({len(in_meta_not_blast)}):")
    for s in sorted(in_meta_not_blast):
        print(f"   {s}")
else:
    print("✅  All metadata sequences found in BLAST.")

# ── Stop here and fix mismatches before proceeding ──────────────────────────
assert not in_blast_not_meta, "Fix missing metadata entries before plotting."
