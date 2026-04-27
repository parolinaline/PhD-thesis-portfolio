getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/BEAST/BEAST_fig")

# BiocManager::install("YuLab-SMU/treedataverse")
#install.packages("ggplot2")
#install.packages("remotes")
#remotes::install_github("YuLab-SMU/ggtree")
#remotes::install_github("YuLab-SMU/treeio")
#install.packages("ggnewscale")
#install.packages("aplot")

#library(treedataverse)
library(ggplot2)
packageVersion("ggplot2")  # should be 4.0.2
library(treeio)
library(ggtree)
library(ggnewscale)
library(aplot)
options("aplot_guides" = "keep")
library(lubridate)
library(magick)
library(ggstance) # For geom_barh
library(ggsci)
library(rhierbaps)


### Following installed with treedataverse
#library(ggtree)
#library(treeio)
#library(ggplot2)
#library(ggtreeExtra)

###############################
#### RHIERBAPS CLUSTERS #######
###############################

snp_matrix <- load_fasta("input/ST654_GIVE_ALIGN.fasta")
results <- hierBAPS(snp_matrix, max.depth = 2, n.pops = 20)
write.csv(results$partition.df, "rhierbaps_clusters.csv", row.names = FALSE)

###########################
#### tree for paper#######
##########################

### Use Figtree and import txt file - add colour using selection mode = Clade
### and annotate nodes from tips command in Tree tab
beast_tree <- read.beast("input/figtree_give_strict.nexus")


### add data - note need to name variables differently to those used in Figtree to label nodes from tips
tip_beast<-data.frame(read.csv("input/metadata_giveST654.csv"))

tip_beast$IncI <- as.character(tip_beast$IncI)

# Merge data with tree
annotated_beast_tree <- full_join(beast_tree, tip_beast, by=c("label"="label"))

### check basic tree
ggtree(annotated_beast_tree, mrsd="2022-11-15", aes(color=Clades)) + 
  theme_tree2()

### Start creating composite figure
p<- ggtree(annotated_beast_tree, mrsd="2022-11-15", aes(color=Clades)) + 
  ylim(-5, 70) + theme_tree2() # +
# scale_x_continuous(sec.axis = dup_axis())

### add vertical lines to demarcate sampling period - could use shading instead
p<- p + geom_vline(xintercept = 2005, linetype = "dotted", color = "darkgrey")
p<- p + geom_vline(xintercept = 2019, linetype = "dotted", color = "darkgrey")
p<- p + geom_vline(xintercept = 2022, linetype = "dotted", color = "darkgrey")

p
### create files for heatmaps
heat<-data.frame(tip_beast$IncI)
heat2<-data.frame(tip_beast$line) ## line through heatmaps
heat3<-data.frame(tip_beast$Source)

### create file for adding vertical lines through heatmaps
line1<-data.frame(tip_beast$full_len)

### link to tree
rownames(heat)<-tip_beast$label
rownames(heat2)<-tip_beast$label
rownames(heat3)<-tip_beast$label

### specify column names
colnames(heat)<-c("IncI")
colnames(heat2)<-c(" ")
colnames(heat3)<-c("Source")

tip_beast$IncI <- as.character(tip_beast$IncI)

### heatmaps

p2 <- gheatmap(p, heat, offset=1.5, width=.02, colnames=T, color=NA, 
               font.size=3, colnames_position="top", colnames_angle=45, 
               colnames_offset_y=1, hjust=0) +
  scale_fill_manual(values=c("0"="white", "1"="black"),
                    name="IncI plasmid",
                    labels=c("0"="Absent", "1"="Present"),
                    na.value="grey90") +
  new_scale_fill()

### vertical lines (no legend needed for this layer)
p2 <- gheatmap(p2, heat2, offset=2.5, width=0.001, colnames=T, color=NA, font.size=0) +
  scale_fill_gradient(low="white", high="black", guide="none") +
  new_scale_fill()

### add Source heatmap
cols <- c(Feline='#ff7d00', Bovine='#8f00ff', Ovine='#ffdd00', Equine='#01befe', Human='#adff02', Canine='#ff006d')
p2 <- gheatmap(p2, heat3, offset=4, width=0.02, colnames=T, color=NA, font.size=3, 
               colnames_position="top", colnames_angle=45, colnames_offset_y=1, hjust=0) +
  scale_fill_manual(values=cols, name="Source")
p2
### Also add a legend for the tree branches (Clades colour)
p3 <- p2 + 
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))
p3


ggsave("output/time_tree_give.png", p3, width=15, height=16, units="cm")



#==============================================================================#
#==========================  PART 02: ZOOMING  ================================#
#==============================================================================#

#### zoom for inset plot

beast_tree <- read.beast("input/figtree_give_strict.nexus")
### add data
tip_beast <- data.frame(read.csv("input/metadata_giveST654.csv"))
# Merge data with tree
annotated_beast_tree <- full_join(beast_tree, tip_beast, by = c("label" = "label"))


p <- ggtree(annotated_beast_tree, mrsd = "2022-11-15", aes(color = Clades)) +
  theme_tree2()

p <- p + xlim(1975, 2023)

p <- p + theme(legend.position = "right")

### Demarcate bounds of sampling period
p<- p + geom_vline(xintercept = 2005, linetype = "dotted", color = "darkgrey")
p<- p + geom_vline(xintercept = 2019, linetype = "dotted", color = "darkgrey")
p<- p + geom_vline(xintercept = 2022, linetype = "dotted", color = "darkgrey")


### Read in time series of case data per clade
ts <- read.csv("input/give_time_series.csv")
ts <- as.data.frame(ts)
ts$A <- as.numeric(ts$A)
ts$B <- as.numeric(ts$B)
ts$C <- as.numeric(ts$C)
ts$D <- as.numeric(ts$D)
ts$Date <- as.Date(ts$Date, "%d/%m/%Y")
ts$Date2 <- decimal_date(ts$Date)


# TODO: update colours to match your tree clade colours
inset_A <- ggplot(ts, aes(x = Date2, y = A, group = 1)) +
  geom_line(colour = "red", linewidth = 0.8) +
  ylab("Clade A") +
  ylim(0, 3) +
  xlim(2005, 2023) +
  ggfun::theme_noxaxis(axis.line = element_line(linetype = 0)) +
  labs(x = NULL)

inset_B <- ggplot(ts, aes(x = Date2, y = B, group = 1)) +
  geom_line(colour = "#739559", linewidth = 0.8) +
  ylab("Clade B") +
  ylim(0, 3) +
  xlim(2005, 2023) +
  ggfun::theme_noxaxis(axis.line = element_line(linetype = 0)) +
  labs(x = NULL)

inset_C <- ggplot(ts, aes(x = Date2, y = C, group = 1)) +
  geom_line(colour = "#39dbff", linewidth = 0.8) +
  ylab("Clade C") +
  ylim(0, 3) +
  xlim(2005, 2023) +
  ggfun::theme_noxaxis(axis.line = element_line(linetype = 0)) +
  labs(x = NULL)

inset_D <- ggplot(ts, aes(x = Date2, y = D, group = 1)) +
  geom_line(colour = "purple", linewidth = 0.8) +
  ylab("Clade D") +
  ylim(0, 3) +
  xlim(2005, 2023) +
  ggfun::theme_noxaxis(axis.line = element_line(linetype = 0)) +
  labs(x = NULL)

### Read in effective population size data
tsNe <- read.csv("input/give_popsize_skyline.csv")  # TODO: update filename
tsNe$Ne  <- as.numeric(tsNe$median)
tsNe$LCI <- as.numeric(tsNe$lower)
tsNe$UCI <- as.numeric(tsNe$upper)


inset_Ne <- ggplot(tsNe, aes(x = time, y = Ne)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "darkgrey", alpha = 0.3) +
  ylab("Ne \n 95% CrI") +
  # ylim(0, 142) +   # TODO: adjust y limit based on your Ne values
  xlim(2005, 2023) +
  labs(x = NULL)



#==============================================================================#
#=======================  PART 03: ADD TIME SERIES  ===========================#
#==============================================================================#

### Add time series plots to tree
p2 <- p %>%
  insert_bottom(inset_A, height = 0.15) %>%
  insert_bottom(inset_B, height = 0.15) %>%
  insert_bottom(inset_C, height = 0.15) %>%
  insert_bottom(inset_D, height = 0.15) %>%
  insert_bottom(inset_Ne, height = 0.3)

p2

ggsave("output/give_time_tree_zoom.png", p2, width = 20, height = 30, units = "cm")

# Define the paths to your images for composite image
base_image_path   <- "output/time_tree_give.png"
overlay_image_path <- "output/give_time_tree_zoom.png"
output_image_path  <- "output/give_time_tree_zoom_overlay.png"

# Read the base image
base_image <- image_read(base_image_path)

# Read the image to be inserted
overlay_image <- image_read(overlay_image_path)

# Resize the overlay image if necessary
overlay_image_resized <- image_scale(overlay_image, "1050x1050")

# Composite the overlay image onto the base image
combined_image <- image_composite(base_image, overlay_image_resized, offset = "+50+50")

# Write the combined image to a new PNG file
image_write(combined_image, path = output_image_path, format = "png")
