getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/BEAST/Give/")



library(ape)
library(lubridate)

tree <- read.tree("iqtree_gubbinsclair_Give.treefile")
tree <- read.tree("gubbins_Give_variantsites.tre")
tree <- read.tree("ST654.treefile")

# Load metadata
meta <- read.csv("Metadata_isolates.csv")

# Convert from DD-MM-YYYY to YYYY-MM-DD
meta$date <- format(dmy(meta$date), "%Y-%m-%d")

# Get ST654 isolates that ALSO exist in the tree
st654_ids <- meta$id[meta$ST == 654]
st654_ids <- intersect(st654_ids, tree$tip.label)

# Prune to keep only ST654
tree_st654 <- keep.tip(tree, st654_ids)

# Create a named vector of dates for easy lookup
dates <- setNames(meta$date, meta$id)

# Rename tips: id|date
tree_st654$tip.label <- sapply(tree_st654$tip.label, function(x) {
  paste0(x, "|", dates[x])
})

write.tree(tree_st654, "gubbins_GiveST654_YMD_29mar.tre")

# Optional: check how many were kept vs filtered
message(length(st654_ids), " ST654 isolates found in tree")

###############################################################################
##############   Do the same but with the alignment file ######################
###############################################################################

# Load the alignment
aln <- read.FASTA("gubbinsclair_Give_polymorphic_sites.fasta")
aln <- read.FASTA("st654_only_gubbins.fasta")

# Get ST654 IDs that exist in BOTH metadata and alignment
st654_ids_aln <- meta$id[meta$ST == 654]
st654_ids_aln <- intersect(st654_ids_aln, names(aln))

# Keep only ST654 sequences that exist in the alignment
aln_st654 <- aln[st654_ids_aln]

# Rename sequences: id|date (same format as tree)
names(aln_st654) <- sapply(names(aln_st654), function(x) {
  paste0(x, "|", dates[x])
})

write.FASTA(aln_st654, "Gubbins_ST654_YMD_29mar.fasta")

message(length(aln_st654), " sequences written")
