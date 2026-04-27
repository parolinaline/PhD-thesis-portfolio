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

# ── 2. Load metadata ─────────────────────────────────────────────────────────
meta = pd.read_excel("./IncFIB_plasmids_information.xlsx")
meta["label"] = meta["FILE"].apply(clean_name)

# Validate all BLAST sequences are in metadata
all_seqs = set(blast["query"]) | set(blast["subject"])
missing  = all_seqs - set(meta["label"])
assert not missing, f"Missing from metadata: {missing}"

# Map isolate label → cluster
label_to_cluster = meta.set_index("label")["Cluster_cdHit"].to_dict()

blast["query_cluster"]   = blast["query"].map(label_to_cluster)
blast["subject_cluster"] = blast["subject"].map(label_to_cluster)

# ── 3. Build cluster-averaged identity matrix ────────────────────────────────
# Keep only inter-isolate comparisons (self-hits already removed upstream)
# Weight by alignment length when averaging within a cluster pair
cluster_grouped = blast.groupby(["query_cluster", "subject_cluster"]).apply(
    lambda g: (g["pident"] * g["length"]).sum() / g["length"].sum()
).reset_index(name="avg_pident")

all_clusters = sorted(set(
    cluster_grouped["query_cluster"].tolist() +
    cluster_grouped["subject_cluster"].tolist()
))
cluster_labels = [f"Cluster {c}" for c in all_clusters]

# Build symmetric matrix; diagonal = 100 (within-cluster identity)
matrix = pd.DataFrame(np.nan, index=cluster_labels, columns=cluster_labels)
np.fill_diagonal(matrix.values, 100.0)

for _, row in cluster_grouped.iterrows():
    q = f"Cluster {int(row['query_cluster'])}"
    s = f"Cluster {int(row['subject_cluster'])}"
    val = row["avg_pident"]
    matrix.loc[q, s] = val
    matrix.loc[s, q] = val

# ── 4. Build binary ST presence matrix ───────────────────────────────────────
# Clean ST values: drop "-" and "nan"
meta_clean = meta[~meta["ST"].astype(str).isin(["-", "nan"])].copy()
meta_clean["ST"] = meta_clean["ST"].astype(str)

all_sts = sorted(meta_clean["ST"].unique(), key=lambda x: int(x) if x.isdigit() else x)
st_labels = [f"ST{st}" for st in all_sts]

# One row per cluster, one column per ST: 1 if present, 0 if not
st_presence = pd.DataFrame(0, index=cluster_labels, columns=st_labels)
for cluster_id in all_clusters:
    label = f"Cluster {cluster_id}"
    sts_in_cluster = meta_clean[
        meta_clean["Cluster_cdHit"] == cluster_id
    ]["ST"].unique()
    for st in sts_in_cluster:
        st_presence.loc[label, f"ST{st}"] = 1

# ── 5. Add n sequences per cluster as annotation ─────────────────────────────
n_seqs = meta.groupby("Cluster_cdHit").size().reset_index(name="n")
n_seqs["label"] = n_seqs["Cluster_cdHit"].apply(lambda x: f"Cluster {x}")
n_seqs = n_seqs.set_index("label")["n"]

# ── 6. Figure layout ─────────────────────────────────────────────────────────
n_clusters = len(cluster_labels)
n_sts      = len(st_labels)

# width ratios: [identity heatmap | ST binary heatmap]
# Identity heatmap is square; ST panel width scales with number of STs
id_ratio = n_clusters
st_ratio = max(2, n_sts * 0.8)

fig, axes = plt.subplots(
    1, 2,
    figsize=(10 + n_sts * 0.4, max(5, n_clusters * 0.55)),
    gridspec_kw={"width_ratios": [id_ratio, st_ratio], "wspace": 0.05}
)
ax_id, ax_st = axes

# ── 7. Identity heatmap ───────────────────────────────────────────────────────
off_diag = matrix.values[~np.eye(n_clusters, dtype=bool)]
vmin = np.nanmin(off_diag) if not np.all(np.isnan(off_diag)) else 0

sns.heatmap(
    matrix.astype(float),
    annot=True, fmt=".1f", annot_kws={"size": 8},
    cmap="YlGnBu",
    vmin=vmin, vmax=100,
    linewidths=0.5, linecolor="white",
    square=True,
    cbar_kws={"label": "Avg % identity (weighted)", "shrink": 0.6, "pad": 0.02},
    ax=ax_id
)
ax_id.set_title("Average pairwise BLASTn identity\nbetween pVir plasmid clusters",
                fontsize=11, fontweight="bold", pad=10)
ax_id.set_xlabel("")
ax_id.set_ylabel("")
ax_id.tick_params(axis="x", rotation=45, labelsize=9)
ax_id.tick_params(axis="y", rotation=0,  labelsize=9)

# Add n sequences as a parenthetical to y-axis tick labels
new_ylabels = [f"{lbl}  (n={n_seqs.get(lbl, '?')})" for lbl in cluster_labels]
ax_id.set_yticklabels(new_ylabels, fontsize=9)

# ── 8. Binary ST presence heatmap ────────────────────────────────────────────
sns.heatmap(
    st_presence.astype(float),
    annot=False,
    cmap=["white", "black"],   # 0 = white, 1 = black
    vmin=0, vmax=1,
    linewidths=0.5, linecolor="#cccccc",
    square=True,
    cbar=False,
    ax=ax_st
)
ax_st.set_title("ST presence\nper cluster",
                fontsize=11, fontweight="bold", pad=10)
ax_st.set_xlabel("")
ax_st.set_ylabel("")
ax_st.tick_params(axis="x", rotation=45, labelsize=9)
ax_st.set_yticks([])   # y labels already on identity heatmap

# ── 9. Save ───────────────────────────────────────────────────────────────────
plt.savefig("./heatmap_pvir_clusters_annotated.png",
            dpi=300, bbox_inches="tight")
plt.savefig("./heatmap_pvir_clusters_annotated.pdf",
            bbox_inches="tight")
print(f"Done! {n_clusters} clusters × {n_clusters} identity matrix, {n_sts} STs in binary panel.")
