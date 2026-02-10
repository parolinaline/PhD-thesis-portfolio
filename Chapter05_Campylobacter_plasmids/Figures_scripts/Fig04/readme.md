#### Description

Figure 4  – Phylogenetic relationships and plasmid marker profiles of Campylobacter jejuni and C. coli isolates from New Zealand. A core genome MLST (cgMLST) tree was constructed to illustrate the phylogenetic relationships among the isolates. Tip colours indicate species (C. jejuni in green and C. coli in pink), and tip shapes represent the host source (bovine, human, poultry, or ovine). Adjacent to the tree, a heatmap displays the sequence type (ST) assigned to each isolate. UD = undetermined, indicating a combination of MLST alleles that has no ST assigned. The coloured tiles of the last heatmap show the presence and type of plasmid markers detected. Each colour indicates the type of plasmid mobilisation predicted by MOB-suite or the complete absence of plasmid sequences as follows: conjugative (pink), mobilisable (purple), non-mobilisable (blue), no plasmid sequences were detected by RFPlasmid (yellow), or absence of that respective plasmid type (grey). ND = Not detected, indicating there were no plasmid markers detected on the plasmid sequence carried by the isolate.

### Input

We are going to use the cgMLST tree obtained as described in methods:  
[cgMLST](../../Methods/cgMLST.md)  

Input files necessary to annotate the tree with the script below:  
Output cgMLST tree:

[Tree](./cgMLST_grape.nwk)
[Relabel the tree](./new_tree_labels.csv)  
[Isolates metadata](./isolates_metadata.csv)
[Plasmid typing results](./bi_heatmap_plasmid.xlsx)  


### Tree annotation

Now let's annotate the tree using R:

```r
getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/Writing_Chapter05/Thesis_Figures_scripts_data/Fig04_cgMLST/")

library(ape)
library(tidytree)
library(ggtree)
library(dplyr)
library(readxl)
library(reshape2)
library(ggplot2)
library(readr)
library(ggfun)
library(aplot)
library(tidyr)
library(treeio)
library(tidyverse)

############################### TREE & LABELS ##################################
# Load tree
phylo_tree <- read.tree("cgMLST_grape.nwk")


###################### RELABING MY TREE   ######################################
###################### EXTRACT THE CURRENT NAME ################################
#write.csv(data.frame(Old = phylo_tree$tip.label), "current_tree_labels.csv", row.names = FALSE)


##################### PREPARE THE CSV FILE WITH THE OLD NAME AND THE NEW NAME ##
# Read the label mapping file ##################################################
label_map <- read.csv("new_tree_labels.csv", stringsAsFactors = FALSE)

# Check structure
head(label_map)


# Replace the tip labels
# Make sure the order of the new labels matches the order in the tree
label_map <- label_map[match(phylo_tree$tip.label, label_map$Old), ]
all(label_map$Old == phylo_tree$tip.label)  # Should be TRUE

# Apply new labels
phylo_tree$tip.label <- label_map$New


################################################################################
############### SPECIES ANNOTATION  ############################################

#Annotation info for the tree
info <- read.csv("./isolates_metadata.csv", header = TRUE)
species_colors <- c("coli" = "#C582B2", "jejuni" = "#67bd57")  # Adjust with your species names and colors

str(info)
all(info$ID == phylo_tree$tip.label)  
info <- info[match(phylo_tree$tip.label, info$ID), ]
all(info$ID == phylo_tree$tip.label) # Should return TRUE

#Find tip labels in the tree that are not in your metadata
setdiff(phylo_tree$tip.label, info$ID)
#Find metadata IDs not present in the tree
setdiff(info$ID, phylo_tree$tip.label)


# Visualize the tree with tip labels
p_tree <- ggtree(phylo_tree) %<+% info +
  geom_tippoint(aes(color = Species, shape = Host), size = 3) +  # Species color and Host shape
#  geom_tiplab(aes(label = label), size = 3, hjust = -0.1) +      # Add tip labels
  scale_color_manual(values = species_colors,
                     labels = c(expression(italic("C. coli")), expression(italic("C. jejuni")))) +
  scale_shape_manual(
    values = c(
      "Bovine" = 8,  # square with +
      "Human" = 17,   # triangle
      "Poultry" = 15, # filled square
      "Ovine" = 19    # filled circle
    ),
    labels = c("Bovine", "Human", "Poultry", "Ovine"),
    name = "Source"
  ) +
  theme_tree() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

p_tree



######################  HEATMAP   ##############################################
#####################  PLASMID MARKERS    ######################################
################################################################################

#HEATMAP WITH THE PLASMID MARKERS COMBINATION PER CONTIG

data_long <- read_excel("bi_heatmap_plasmid.xlsx", sheet = "long_format")


# 3️⃣ Define plasmid groups
conjugative <- c("M04_R08_R37_R59_S03_V01","M04_R59_S03_V01","M05_R02_S01_V03",
                 "M05_R05_S01_V03","M10_R10_S02_V02","R59_S03_V01")
mobilisable <- c("M06","M06_M06_R03_R04","M06_R03","M17_R04")
non_mobilisable <- c("ND","R04","R08","R10","R11","R15","R36","R59")
no_plasmids <- c("No_plasmids")

# 4️⃣ Add mobility type & color
data_long <- data_long %>%
  mutate(
    Mobility_type = case_when(
      Presence == 1 & Plasmid_marker %in% conjugative ~ "Conjugative",
      Presence == 1 & Plasmid_marker %in% mobilisable ~ "Mobilisable",
      Presence == 1 & Plasmid_marker %in% non_mobilisable ~ "Non-mobilisable",
      Presence == 1 & Plasmid_marker %in% no_plasmids ~ "No_plasmids",
      TRUE ~ "Absence"
    ),
    Color = case_when(
      Mobility_type == "Conjugative" ~ "#DE3163",
      Mobility_type == "Mobilisable" ~ "purple",
      Mobility_type == "Non-mobilisable" ~ "steelblue",
      Mobility_type == "No_plasmids" ~ "#BFB040",
      Mobility_type == "Absence" ~ "#E5E4E2"
    )
  )




# 5️⃣ Plot with legend
heatmap_plot <- ggplot(data_long, aes(x = Plasmid_marker, y = ID, fill = Mobility_type)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "Conjugative" = "#DE3163",
      "Mobilisable" = "purple",
      "Non-mobilisable" = "steelblue",
      "No_plasmids" = "#BFB040",
      "Absence" = "#E5E4E2"
    ),
    name = "Plasmid mobility"  # legend title
  ) +
  labs(title = "Plasmid types") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

heatmap_plot



######################### ST HEATMAP ###########################################
################################################################################
# Convert ST to factor (helps keep consistent order and colors)
library(ggplot2)
library(randomcoloR)

# Convert ST to factor (helps keep consistent order and colors)
info$ST <- as.factor(info$ST)

# Automatically generate a distinct color for each ST
num_STs <- length(unique(info$ST))
st_colors <- distinctColorPalette(num_STs)
names(st_colors) <- levels(info$ST)

# Build the one-column heatmap with ST labels inside tiles
st_heatmap <- ggplot(info, aes(x = "ST", y = ID, fill = ST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ST), size = 3, color = "black") +  # add text labels
  scale_fill_manual(values = st_colors) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"  # remove legend
  )

st_heatmap

##########################BAR PLOT##############################################
########################PLASMID SIZE############################################
################################################################################
plasmid_size <- read.csv("barplot.csv")

plasmid_size

gbar <- plasmid_size  # already in a suitable format
# Create the bar plot
pbar <- ggplot(gbar, aes(x = ID, y = Length)) +
  geom_bar(stat = "identity", fill = "#1F6683")+
  xlab(NULL) +
  coord_flip() +
  ggfun::theme_noyaxis() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )


pbar


######################### APLOT ################################################
################# PLOTTING EVERYTHING TOGETHER #################################
################################################################################

final_plot_contigs <- pbar %>%
  insert_left (heatmap_plot, width = 0.7) %>%
  insert_left(st_heatmap, width = 0.1) %>%
  insert_left(p_tree, width = 0.8)

final_plot_contigs

#save like 10 x 14

```
