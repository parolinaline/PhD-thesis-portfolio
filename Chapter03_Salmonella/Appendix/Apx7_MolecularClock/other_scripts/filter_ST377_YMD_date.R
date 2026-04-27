getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/BEAST_prep/bovis/")
getwd()


library(ape)
library(lubridate)

tree <- read.tree("iqtree_gubbins_masked_bovis.treefile")
tree <- read.tree("gubbins_Bovis_variantsites.tre")

# Load metadata
meta <- read.csv("Metadata_isolates.csv")

# Convert from DD-MM-YYYY to YYYY-MM-DD
meta$date <- format(dmy(meta$date), "%Y-%m-%d")

# Get ST377 isolates that ALSO exist in the tree
st377_ids <- meta$id[meta$ST == 377]
st377_ids <- intersect(st377_ids, tree$tip.label)

# Prune to keep only ST377
tree_st377 <- keep.tip(tree, st377_ids)

# Create a named vector of dates for easy lookup
dates <- setNames(meta$date, meta$id)

# Rename tips: id|date
tree_st377$tip.label <- sapply(tree_st377$tip.label, function(x) {
  paste0(x, "|", dates[x])
})

write.tree(tree_st377, "gubbins_ST377bovis_YMD_variantsites.tre")

# Optional: check how many were kept vs filtered
message(length(st377_ids), " ST377 isolates found in tree")

###############################################################################
##############   Do the same but with the alignment file ######################
###############################################################################

# Load the alignment
aln <- read.FASTA("gubbins_bovis.filtered_polymorphic_sites.fasta")

# Get ST377 IDs that exist in BOTH metadata and alignment
st377_ids_aln <- meta$id[meta$ST == 377]
st377_ids_aln <- intersect(st377_ids_aln, names(aln))

# Keep only ST377 sequences that exist in the alignment
aln_st377 <- aln[st377_ids_aln]

# Rename sequences: id|date (same format as tree)
names(aln_st377) <- sapply(names(aln_st377), function(x) {
  paste0(x, "|", dates[x])
})

write.FASTA(aln_st377, "Gubbins_ST377_YMD_align.fasta")

message(length(aln_st377), " sequences written")
