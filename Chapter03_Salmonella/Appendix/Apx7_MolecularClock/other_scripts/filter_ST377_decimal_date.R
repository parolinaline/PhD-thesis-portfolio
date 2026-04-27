getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/BEAST_prep/bovis/")
getwd()

library(ape)
library(lubridate)

tree <- read.tree("iqtree_gubbins_masked_bovis.treefile")

# Load metadata
meta <- read.csv("Metadata_isolates.csv")

# Convert dates to decimal format for BEAST
# Assuming your dates are in DD-MM-YYYY or similar format
meta$decimal_date <- decimal_date(dmy(meta$date))

# Get ST377 isolates that ALSO exist in the tree
st377_ids <- meta$id[meta$ST == 377]
st377_ids <- intersect(st377_ids, tree$tip.label)

# Prune to keep only ST377
tree_st377 <- keep.tip(tree, st377_ids)

# Create a named vector of DECIMAL dates for easy lookup
dates <- setNames(meta$decimal_date, meta$id)

# Rename tips: id|decimal_date
tree_st377$tip.label <- sapply(tree_st377$tip.label, function(x) {
  paste0(x, "|", dates[x])
})

write.tree(tree_st377, "iqtreeclair_ST377_decimalY.tre")

# Optional: check how many were kept vs filtered
message(length(st377_ids), " ST377 isolates found in tree")

###############################################################################
##############   Do the same but with the alignment file ######################
###############################################################################

# Load the Gubbins masked alignment
aln <- read.FASTA("gubbins_bovis.filtered_polymorphic_sites.fasta")

# Keep only ST377 sequences
aln_st377 <- aln[names(aln) %in% st377_ids]

# Rename headers: id|decimal_date (dates vector already has decimal dates)
names(aln_st377) <- sapply(names(aln_st377), function(x) {
  paste0(x, "|", dates[x])
})

write.FASTA(aln_st377, "ST377_decimalyear_gubbinsclair.fasta")
message(length(aln_st377), " sequences written")


