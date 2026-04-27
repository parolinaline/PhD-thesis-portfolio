
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/4_Fig_Give_plasmids/")

getwd()
# Load required libraries
library(ape)
library(tidytree)
library(ggtree)
library(dplyr)
library(readxl)
library(reshape2)
library(ggplot2)
library(ggfun)
library(aplot)
library(tidyr)
library(treeio)
library(tidyverse)
library(ggtext)
library(extrafont)
library(viridis)
library(ggnewscale)
loadfonts(device = "win")


############################### TREE & LABELS ##################################
# Load tree
phylo_tree <- read.tree("Give_plasmids_roary_tree.newick")

###################### RELABING MY TREE
###################### EXTRACT THE CURRENT NAME
#write.csv(data.frame(Old = phylo_tree$tip.label), "current_tree_labels.csv", row.names = FALSE)

##################### PREPARE THE CSV FILE WITH THE OLD NAME AND THE NEW NAME
# Read the label mapping file
label_map <- read.csv("new_tree_labels_give.csv", stringsAsFactors = FALSE)

# Check structure
head(label_map)

# Replace the tip labels
# Make sure the order of the new labels matches the order in the tree
label_map <- label_map[match(phylo_tree$tip.label, label_map$Old), ]
all(label_map$Old == phylo_tree$tip.label)  # Should be TRUE

# Apply new labels
phylo_tree$tip.label <- label_map$New

#=============================================================================#
############################## METADATA + TREE ################################
#=============================================================================#

# Load metadata
info <- read.csv("./roary_annotation_metadata_Give.csv", header = TRUE)
#info <- info[match(phylo_tree$tip.label, info$ID), ]
#stopifnot(all(info$ID == phylo_tree$tip.label))

# Source colours
source_colours <- c("Bovine" = "#A65DD0", "Human" = "#45BA81", "Poultry" = "#5C62D1", 
                    "Canine" = "#F0FC03", "Equine" = "#B51C1C")

str(info)
all(info$ID == phylo_tree$tip.label)  # Should return TRUE
info <- info[match(phylo_tree$tip.label, info$ID), ]
all(info$ID == phylo_tree$tip.label)

##### 1) after reading tree & relabeling tips #####
# keep tip order in a vector
tip_levels <- phylo_tree$tip.label

##### 2) ensure metadata matches tree BEFORE plotting #####
# (you already do this but keep the check)
info <- info[match(tip_levels, info$ID), ]
stopifnot(all(info$ID == tip_levels))  # will stop early if mismatch

##### 3) tidy up p_Give #####
p_Give <- ggtree(phylo_tree) %<+% info + 
  geom_tippoint(aes(color = Source, shape = Mob_prediction), size = 3.8) +
#  geom_tiplab(size = 3, offset = 0.5, hjust = 0) +
  scale_color_manual(values = source_colours, labels = c("Bovine","Human","Poultry","Canine","Equine")) +
  scale_shape_manual(
    values = c("conjugative" = 16, "mobilisable" = 17, "non-mobilisable" = 15),
    labels = c("Conjugative", "Mobilisable", "Non-mobilisable"),
    name = "Mobility prediction"
  ) +
  theme_tree() +
  theme(legend.text = element_text(size=12), legend.title = element_text(size = 12))

p_Give

#=============================================================================#
########################### PLASMID MARKER HEATMAP ###########################
#=============================================================================#

# Example: plasmid markers heatmap
data_pm <- read.csv("binary_heatmap_plasmids_Give.csv", row.names = 1, check.names = FALSE)
data_long_pm <- reshape2::melt(as.matrix(data_pm))
colnames(data_long_pm) <- c("Isolate", "Plasmid_markers", "Presence")

# IMPORTANT: set factor levels (try tip_levels or rev(tip_levels) if rows are flipped)
data_long_pm$Isolate <- factor(data_long_pm$Isolate, levels = tip_levels)

nanopore_heatmap <- ggplot(data_long_pm, aes(x = Plasmid_markers, y = Isolate, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#D04E59"), name = "Plasmid types", labels= c("Absent", "Present")) +
  theme_minimal() +
#  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab(NULL)

nanopore_heatmap


#=============================================================================#
########################### PLASMID MARKER AMR GENES ##########################
#=============================================================================#

# AMR heatmap (same idea)
data_amr <- read.csv("binary_heatmap_AMR_genes.csv", row.names = 1, check.names = FALSE)
data_long_amr <- reshape2::melt(as.matrix(data_amr))
colnames(data_long_amr) <- c("Isolate", "AMR_Genes", "Presence")
data_long_amr$Isolate <- factor(data_long_amr$Isolate, levels = tip_levels)

AMR_heatmap <- ggplot(data_long_amr, aes(x = AMR_Genes, y = Isolate, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#088F8F"), name = "AMR markers", labels = c("Absent", "Present")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.y = element_blank()) +
  ylab(NULL)

AMR_heatmap


#=============================================================================#
#################### ONE COLUMN HEATMAP BANDAGE RESULTS #######################
#=============================================================================#

# bandage heatmap (ensure ID factor)
data_bandage <- read.csv("bandage_results.csv")
data_bandage$Circularisable <- as.factor(data_bandage$Circularisable)
data_bandage$ID <- factor(data_bandage$ID, levels = tip_levels)

nanopore_heatmap_bandage <- ggplot(data_bandage, aes(x = "Circularisable", y = ID, fill = Circularisable)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "black"), name = "Circularisable", labels = c("No","Yes")) +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


nanopore_heatmap_bandage


#=============================================================================#
########################### PLASMID SIZE BAR PLOT #############################
#=============================================================================#

# bar plot: set ID factor with same levels (and ensure same orientation as tree)
plasmid_size <- readxl::read_xlsx("plasmid_size_annotation.xlsx", sheet = "plasmid_size")
plasmid_size$ID <- factor(plasmid_size$ID, levels = tip_levels)

pbar_size <- ggplot(plasmid_size, aes(x = ID, y = Length)) +
  geom_bar(stat = "identity", fill = "#1F6683") +
  coord_flip() +
  theme_minimal() +
  ggfun::theme_noyaxis() +
  xlab(NULL) +
  ylab(NULL)

pbar_size


#=============================================================================#
############################## DUMBBELL GC PLOT  ##############################
#=============================================================================#

# GC (dumbbell) plot: ensure ID factor levels too
gc_info <- read.csv("./gc_content.csv", header = TRUE)
gc_info$ID <- factor(gc_info$ID, levels = tip_levels)

data_gc <- gc_info%>%
  select(ID, gc_chromosome, gc_plasmid) %>%
  mutate(diff = gc_chromosome - gc_plasmid) %>%
  pivot_longer(cols = c(gc_chromosome, gc_plasmid)) %>%
  rename( Type = name,
          GC_content = value)
head(data_gc)

#In order to plot the range between the two groups later on, we need two tibbles which are only holding the data of 1 data (gc content from plasmid and from chromosome) each.
plasmid <- data_gc %>%
  filter( Type == "gc_plasmid")
head(plasmid)

chromosome <- data_gc %>%
  filter( Type == "gc_chromosome")
head(chromosome)


# let's calculate the mean from the GC content of plasmids and chromosomes by species

stats_plasmid <- data_gc %>%
  filter(Type == "gc_plasmid") %>%
  summarise(mean = mean(GC_content)) %>%
  pull(mean)  # Extract the value
head(stats_plasmid)

# Calculate mean GC content for chromosome by species
stats_chromosome <- data_gc %>%
  filter(Type == "gc_chromosome") %>%
  summarise(mean = mean(GC_content)) %>%
  pull(mean)  # Extract the value
head(stats_chromosome)

# Now let's put everything together and add the lines showing the mean of the gc content from plasmids and from chromosome
p_gc_content <- ggplot(data_gc)+
  # Add segments first (so they are in the background)
  geom_segment(data = plasmid,
               aes(x = GC_content, y = ID,
                   yend = chromosome$ID, xend = chromosome$GC_content), 
               color = "#E5E4E2",
               size = 4.5, # Note that I sized the segment to fit the points
               alpha = .5) +
  # Add points on top
  geom_point(aes(x = GC_content, y = ID, color = Type), size = 3, show.legend = TRUE)+
  # Color points
  scale_color_manual(values = c("#009688", "#762a83"), name = "GC content type", labels = c("Chromosome", "Plasmid")) +
  # Add theme
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "#4a4e4d"),
        strip.text.y.left  = element_text(angle = 0),
        panel.background = element_rect(fill = "white", color = "white"),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "#4a4e4d", family = "Segoe UI"),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(1,1,.5,1, "cm")) +
  # Add mean GC content lines
  geom_vline(xintercept = stats_plasmid, linetype = "solid", size = 1, alpha = .8, color = "#762a83") +  # Plasmid mean
  geom_vline(xintercept = stats_chromosome, linetype = "solid", size = 1, alpha = .8, color = "#009688")  # chromosome mean

p_gc_content

warnings()


# build your p_final with that gc_info (I omitted full recalc for brevity)

##### 5) Compose starting from the tree and insert_right() #####

final_plot_contigs <- nanopore_heatmap %>%
  insert_left(nanopore_heatmap_bandage, width = 0.08) %>%
  insert_left(p_Give, width = 0.8) %>%
  insert_right(AMR_heatmap, width = 0.6) %>%
  insert_right(pbar_size, width = 0.4) %>%
  insert_right(p_gc_content, width = 0.8)

final_plot_contigs


# now save
ggsave("Give_roarytree_plasmids.svg", plot = final_plot_contigs, width = 22, height = 14, units = "in", dpi = 600)
