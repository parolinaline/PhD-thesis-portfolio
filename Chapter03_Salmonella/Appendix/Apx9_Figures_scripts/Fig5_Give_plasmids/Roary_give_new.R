
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

###################### RELABELLING THE TREE ####################################

# Read the label mapping file
label_map <- read.csv("new_tree_labels_give.csv", stringsAsFactors = FALSE)

# Check structure
head(label_map)

# Replace the tip labels — ensure order matches the tree
label_map <- label_map[match(phylo_tree$tip.label, label_map$Old), ]
all(label_map$Old == phylo_tree$tip.label)  # Should be TRUE

# Apply new labels
phylo_tree$tip.label <- label_map$New

# Save the relabelled tree
write.tree(phylo_tree, file = "Give_NEW.newick")


############################### METADATA #######################################

# Load metadata
info <- read.csv("./roary_annotation_metadata_Give.csv", header = TRUE)

# Source colours
source_colours <- c("Bovine" = "#A65DD0", "Human" = "#45BA81", "Poultry" = "#5C62D1",
                    "Canine" = "#F0FC03", "Equine" = "#B51C1C")

# Match metadata to tree tip order
info <- info[match(phylo_tree$tip.label, info$ID), ]
stopifnot(all(info$ID == phylo_tree$tip.label))

############################### ANNOTATING THE TREE ############################

p_Give <- ggtree(phylo_tree) %<+% info +
  geom_tippoint(aes(color = Source, shape = Mob_prediction), size = 3.8) +
  geom_tiplab(size = 3, offset = 0.5, hjust = 0) +
  scale_color_manual(values = source_colours) +
  scale_shape_manual(
    values = c("conjugative" = 16,
               "mobilisable" = 17,
               "non-mobilisable" = 15),
    labels = c("Conjugative", "Mobilisable", "Non-mobilisable"),
    name = "Mobility prediction"
  ) +
  theme_tree() +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

p_Give

################################################################################
# CRITICAL: Extract tip order from the tree to enforce alignment across panels #
################################################################################

tip_order <- get_taxa_name(p_Give)

############################### HEATMAP ########################################
######################### PLASMID MARKERS ######################################

data <- read.csv("binary_heatmap_plasmids_Give.csv", row.names = 1, check.names = FALSE)

data_long <- melt(as.matrix(data))
colnames(data_long) <- c("ID", "Plasmid_markers", "Presence")

# Enforce tree tip order
data_long$ID <- factor(data_long$ID, levels = tip_order)

# Check for mismatches
cat("Plasmid heatmap — IDs in tree but not in data:", 
    paste(setdiff(tip_order, unique(as.character(data_long$ID))), collapse = ", "), "\n")
cat("Plasmid heatmap — IDs in data but not in tree:", 
    paste(setdiff(unique(as.character(data_long$ID)), tip_order), collapse = ", "), "\n")

nanopore_heatmap <- ggplot(data_long, aes(x = Plasmid_markers, y = ID, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#D04E59"),
                    name = "Plasmid types", labels = c("Absent", "Present")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ylab(NULL)

nanopore_heatmap

############################### HEATMAP ########################################
######################### AMR GENES ############################################

data_amr <- read.csv("binary_heatmap_AMR_genes.csv", row.names = 1, check.names = FALSE)

data_amr_long <- melt(as.matrix(data_amr))
colnames(data_amr_long) <- c("Isolate", "AMR_Genes", "Presence")

# Enforce tree tip order
data_amr_long$Isolate <- factor(data_amr_long$Isolate, levels = tip_order)

# Check for mismatches
cat("AMR heatmap — IDs in tree but not in data:", 
    paste(setdiff(tip_order, unique(as.character(data_amr_long$Isolate))), collapse = ", "), "\n")
cat("AMR heatmap — IDs in data but not in tree:", 
    paste(setdiff(unique(as.character(data_amr_long$Isolate)), tip_order), collapse = ", "), "\n")

AMR_heatmap <- ggplot(data_amr_long, aes(x = AMR_Genes, y = Isolate, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#088F8F"),
                    name = "AMR markers", labels = c("Absent", "Present")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  ylab(NULL)

AMR_heatmap

############################### BAR PLOT #######################################
######################### PLASMID SIZE #########################################

plasmid_size <- read_xlsx("plasmid_size_annotation.xlsx", sheet = "plasmid_size")

gbar <- plasmid_size

# Enforce tree tip order
gbar$ID <- factor(gbar$ID, levels = tip_order)

# Check for mismatches
cat("Bar plot — IDs in tree but not in data:", 
    paste(setdiff(tip_order, unique(as.character(gbar$ID))), collapse = ", "), "\n")
cat("Bar plot — IDs in data but not in tree:", 
    paste(setdiff(unique(as.character(gbar$ID)), tip_order), collapse = ", "), "\n")

pbar <- ggplot(gbar, aes(x = ID, y = Length)) +
  geom_bar(stat = "identity", fill = "#1F6683") +
  xlab(NULL) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

pbar

############################### DUMBBELL PLOT ###################################
######################### GC CONTENT ###########################################

info_gc <- read.csv("./gc_content.csv", header = TRUE)

data_gc <- info_gc %>%
  select(ID, gc_chromosome, gc_plasmid) %>%
  mutate(diff = gc_chromosome - gc_plasmid) %>%
  pivot_longer(cols = c(gc_chromosome, gc_plasmid)) %>%
  rename(Type = name,
         GC_content = value)

# Enforce tree tip order
data_gc$ID <- factor(data_gc$ID, levels = tip_order)

# Separate into plasmid and chromosome subsets (after factoring)
plasmid <- data_gc %>% filter(Type == "gc_plasmid")
chromosome <- data_gc %>% filter(Type == "gc_chromosome")

# Check for mismatches
cat("GC dumbbell — IDs in tree but not in data:", 
    paste(setdiff(tip_order, unique(as.character(data_gc$ID))), collapse = ", "), "\n")
cat("GC dumbbell — IDs in data but not in tree:", 
    paste(setdiff(unique(as.character(data_gc$ID)), tip_order), collapse = ", "), "\n")

# Calculate mean GC content for the vertical reference lines
stats_plasmid <- data_gc %>%
  filter(Type == "gc_plasmid") %>%
  summarise(mean = mean(GC_content, na.rm = TRUE)) %>%
  pull(mean)

stats_chromosome <- data_gc %>%
  filter(Type == "gc_chromosome") %>%
  summarise(mean = mean(GC_content, na.rm = TRUE)) %>%
  pull(mean)

# Build the dumbbell plot
# NOTE: joining plasmid + chromosome to get matched segments avoids row-order issues
segment_data <- inner_join(
  plasmid %>% select(ID, GC_plasmid = GC_content),
  chromosome %>% select(ID, GC_chromosome = GC_content),
  by = "ID"
)

p_final <- ggplot(data_gc) +
  # Segments connecting chromosome and plasmid GC per isolate
  geom_segment(data = segment_data,
               aes(x = GC_plasmid, xend = GC_chromosome, y = ID, yend = ID),
               color = "#E5E4E2", size = 4.5, alpha = 0.5) +
  # Points
  geom_point(aes(x = GC_content, y = ID, color = Type), size = 3, show.legend = TRUE) +
  scale_color_manual(values = c("gc_chromosome" = "#009688", "gc_plasmid" = "#762a83"),
                     name = "GC content type",
                     labels = c("Chromosome", "Plasmid")) +
  # Mean GC reference lines
  geom_vline(xintercept = stats_plasmid, linetype = "solid", size = 1, alpha = 0.8, color = "#762a83") +
  geom_vline(xintercept = stats_chromosome, linetype = "solid", size = 1, alpha = 0.8, color = "#009688") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "#4a4e4d"),
        strip.text.y.left = element_text(angle = 0),
        panel.background = element_rect(fill = "white", color = "white"),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "#4a4e4d"),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(1, 1, 0.5, 1, "cm"))

p_final

############################### HEATMAP ########################################
######################### BANDAGE RESULTS ######################################

data_bandage <- read.csv("bandage_results.csv")

# Enforce tree tip order
data_bandage$ID <- factor(data_bandage$ID, levels = tip_order)
data_bandage$Circularisable <- as.factor(data_bandage$Circularisable)

# Check for mismatches
cat("Bandage heatmap — IDs in tree but not in data:", 
    paste(setdiff(tip_order, unique(as.character(data_bandage$ID))), collapse = ", "), "\n")
cat("Bandage heatmap — IDs in data but not in tree:", 
    paste(setdiff(unique(as.character(data_bandage$ID)), tip_order), collapse = ", "), "\n")

nanopore_heatmap_bandage <- ggplot(data_bandage, aes(x = "Circularisable", y = ID, fill = Circularisable)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "black"),
                    name = "Circularisable", labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank())

nanopore_heatmap_bandage

######################### APLOT ################################################
################# PLOTTING EVERYTHING TOGETHER #################################

# Assemble using the tree as anchor, inserting panels to the right
final_plot_contigs <- p_Give %>%
  insert_right(nanopore_heatmap_bandage, width = 0.08) %>%
  insert_right(nanopore_heatmap, width = 0.8) %>%
  insert_right(AMR_heatmap, width = 0.6) %>%
  insert_right(pbar, width = 0.4) %>%
  insert_right(p_final, width = 0.6)

final_plot_contigs

# Optional: save to file
# ggsave("Give_plasmids_composite_figure.pdf", final_plot_contigs, width = 24, height = 14)
# ggsave("Give_plasmids_composite_figure.png", final_plot_contigs, width = 24, height = 14, dpi = 300)