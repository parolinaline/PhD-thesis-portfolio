
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/5B_Fig_Bovis_plasmids")

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
library(scales)
loadfonts(device = "win")

############################### TREE & LABELS ##################################
# Load tree
phylo_tree <- read.tree("roary_plasmid_bovis.newick")

###################### RELABELING MY TREE
###################### EXTRACT THE CURRENT NAME
write.csv(data.frame(Old = phylo_tree$tip.label), "current_tree_labels_rep08jul.csv", row.names = FALSE)
##################### PREPARE THE CSV FILE WITH THE OLD NAME AND THE NEW NAME

# Read the label mapping file
label_map <- read.csv("new_tree_labels_rep_08july.csv", stringsAsFactors = FALSE)

# Check structure
head(label_map)

# Replace the tip labels
# Make sure the order of the new labels matches the order in the tree
label_map <- label_map[match(phylo_tree$tip.label, label_map$Old), ]
all(label_map$Old == phylo_tree$tip.label)  # Should be TRUE

# Apply new labels
phylo_tree$tip.label <- label_map$New

phylo_tree


############################### HEATMAP DATA ###################################
############################### REPLICON TYPES #################################

info <- read.csv("./Roary_rep_INC.csv", header = TRUE)
info <- info[match(phylo_tree$tip.label, info$ID), ]
stopifnot(all(info$ID == phylo_tree$tip.label))

# Incompatibility Groups or replicon types
INC <- read.csv("./Roary_rep_INC.csv", row.names = 1)
#Stores the row names of heatmapData in a variable called rn. The row names are the first column of the CSV file.
rn <- rownames(INC)
#Converts the values in heatmapData to a character format, so that they can be used as labels for the heatmap.
INC <- as.data.frame(sapply(INC, as.character))
#Assigns the row names that were previously stored in rn as the new row names of heatmapData.
rownames(INC) <- rn

heatmap.colours <- c("1" = "#FFFDD0",
                     "2" = "#adb2fb",
                     "3" = "#a8e6cf",
                     "4" = "red",
                     "5" = "steelblue",
                     "6" = "navy",
                     "7" = "pink",
                     "8" = "magenta")



############################### TREE PLOTTING ##################################

# Base tree with point annotation

p <- ggtree(phylo_tree) +
  theme_tree()+
#  geom_tiplab(size = 3, offset = 0, hjust = 0) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size = 12))
p


################################################################################
########################## Adding Inc one column heatmap ##########################################

p_tree_final <- gheatmap(p, INC, offset = 0, width=0.05, color=NULL, colnames_position = "bottom", 
                         colnames_angle=90, colnames_offset_y = 4, colnames_offset_x = -0.1, 
                         hjust=0.5, font.size=3,legend_title = "Incompatibility Groups") +
  scale_fill_manual(values=heatmap.colours,
                    labels=c("IncFIB", "IncHI2A", "IncI", "IncK2/Z", "IncN", "IncX1", "IncX4", "rep_1760"),
                    name="Plasmid types", breaks=1:8)

p_tree_final


##########################HEATMAP###############################################
#####################   SEQUENCE TYPES      ####################################
################################################################################

#HEATMAP WITH THE ST present within cd-hit clusters
data <- read.csv("binary_heatmap_ST_Bovisplasmids.csv", row.names = 1)

data_long <- melt(as.matrix(data))
colnames(data_long) <- c("Isolate", "Sequence_Type", "Presence")
nanopore_heatmap<-ggplot(data_long, aes(x = Sequence_Type, y = Isolate, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#D04E59"), name = "ST", labels = c("Absent", "Present")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab(NULL)

nanopore_heatmap



##########################HEATMAP###############################################
######################## AMR MARKERS #########################################
################################################################################

#HEATMAP WITH THE PLASMID MARKERS COMBINATION PER CONTIG
AMR_data <- read.csv("binary_heatmap_AMR_Bovisplasmids.csv", row.names = 1)

AMR_data_long <- melt(as.matrix(AMR_data))
colnames(AMR_data_long) <- c("Isolate", "AMR_markers", "Presence")
nanopore_heatmap_AMR<-ggplot(AMR_data_long, aes(x = AMR_markers, y = Isolate, fill = as.factor(Presence))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#4b89bb"), name = "AMR markers", labels = c("Absent", "Present")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.y = element_blank()) +
  ylab(NULL)

nanopore_heatmap_AMR

####################### BOXPLOT #########################################
##################### PLASMID SIZE  ########################################
################################################################################
#Load data
boxplot <- read.csv("boxplot_bovis_rep.csv")


boxplot_figure <- ggplot(boxplot, aes(y = ID, x = Length)) +
  geom_boxplot(fill = "lightblue", color = "darkgrey") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal() +
  labs(x = "Distribution of Plasmid Lengths per Cluster") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

boxplot_figure


###################      BAR PLOT            ###################################
###############       NUMBER OF ISOLATES PER CLUSTER       #####################
################################################################################

no_isolates <- read.csv("number_isolates_per_cluster.csv")
gbar <- no_isolates  # already in a suitable format


pheat <- ggplot(gbar, aes(x = "Isolates", y = ID, fill = Number)) +
  geom_tile(color = "white") +
  # log scale for fill so differences show up
  scale_fill_gradientn(
    trans = "log10",
    colours = c("#D6EAF8", "pink"),
#    breaks = c(1, 20, 100, 300),
    labels = label_number()
  ) +
  # add raw counts inside each tile
  geom_text(aes(label = Number), color = "black", size = 3) +
  theme_minimal() +
  labs(x = NULL, y = NULL, fill = "Number of isolates") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank()
  ) +
  labs(x = "No of isolates per cluster")

pheat


################ADDING A DUMBBELL PLOT WITH THE GC CONTENT######################
##########################DUMBBELL PLOT#########################################
################################################################################

info <- read.csv("./GC_content_BovisPlasmids.csv")

data_gc <- info%>%
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
gc_content_plot <- ggplot(data_gc)+
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
  scale_color_manual(values = c("#009688", "#762a83"),
                     name = "GC content type",
                     labels = c("Chromosome", "Plasmid")) +
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

gc_content_plot

warnings()


######################### APLOT ################################################
################# PLOTTING EVERYTHING TOGETHER #################################
################################################################################
#final_plot_contigs


final_plot_contigs <- nanopore_heatmap_AMR %>% 
  insert_right(boxplot_figure, width = 0.6) %>%
  insert_right(gc_content_plot, width = 0.5) %>%
  insert_left(nanopore_heatmap, width = 0.3) %>% 
  insert_left(pheat, width = 0.1) %>%
  insert_left(p_tree_final, width = 0.8) 
  

final_plot_contigs

################################################################################
################################################################################
# Save the figure as SVG
#install.packages("svglite")
library(svglite)

ggsave("Bovis_rep_plasmids.svg", plot = final_plot_contigs, width = 24, height = 10, units = "in", dpi = 600)


