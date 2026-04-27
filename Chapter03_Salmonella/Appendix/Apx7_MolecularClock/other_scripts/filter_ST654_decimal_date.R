
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/BEAST/Give/")

library(ape)
library(lubridate)

tree <- read.tree("iqtree_gubbinsclair_Give.treefile")
tree <- read.tree("ST654.treefile")
meta <- read.csv("Metadata_isolates.csv")

# Check what's in the tree
message("Total tips in tree: ", length(tree$tip.label))
message("Internal nodes labeled: ", length(tree$node.label))

# Check tip labels match metadata
message("Tips matching metadata: ", sum(tree$tip.label %in% meta$id))
message("Tips NOT in metadata: ", sum(!tree$tip.label %in% meta$id))

# Show any mismatches
if(any(!tree$tip.label %in% meta$id)) {
  message("Missing from metadata: ", paste(tree$tip.label[!tree$tip.label %in% meta$id], collapse=", "))
}

# Convert dates to decimal format
meta$decimal_date <- decimal_date(dmy(meta$date))

# Filter to ST654
st654_ids <- meta$id[meta$ST == 654]
st654_in_tree <- intersect(st654_ids, tree$tip.label)

message("ST654 isolates in metadata: ", length(st654_ids))
message("ST654 isolates found in tree: ", length(st654_in_tree))

# Prune tree
tree_st654 <- keep.tip(tree, st654_in_tree)

# Remove internal node labels (optional, cleaner for BEAST)
tree_st654$node.label <- NULL

# Create named vector of dates
dates <- setNames(meta$decimal_date, meta$id)

# Rename tips with decimal date
tree_st654$tip.label <- sapply(tree_st654$tip.label, function(x) {
  paste0(x, "|", round(dates[x], 4))
})

write.tree(tree_st654, "iqtree_GiveST654_29MAR.tre")

message("Final tree has ", length(tree_st654$tip.label), " tips")
###############################################################################
##############   Do the same but with the alignment file ######################
###############################################################################

# Load the Gubbins masked alignment
aln <- read.FASTA("gubbinsclair_Give_polymorphic_sites.fasta")
aln <- read.FASTA("st654_only_gubbins.fasta")

# Keep only ST654 sequences
aln_st654 <- aln[names(aln) %in% st654_ids]

# Rename headers: id|decimal_date (dates vector already has decimal dates)
names(aln_st654) <- sapply(names(aln_st654), function(x) {
  paste0(x, "|", dates[x])
})

write.FASTA(aln_st654, "ST654give_decimal_29mar.fasta")
message(length(aln_st654), " sequences written")


