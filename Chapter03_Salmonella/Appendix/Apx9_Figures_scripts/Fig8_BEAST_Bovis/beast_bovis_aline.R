getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/BEAST/BEAST_fig")

library(ggplot2)
packageVersion("ggplot2")
library(treeio)
library(ggtree)
library(ggnewscale)
library(aplot)
options("aplot_guides" = "keep")
library(lubridate)
library(magick)
library(ggstance)
library(ggsci)
library(rhierbaps)


###############################
#### RHIERBAPS CLUSTERS #######
###############################

snp_matrix <- load_fasta("input/bovis/Bovis_ST377_decimal_fullalign.fasta")
results <- hierBAPS(snp_matrix, max.depth = 2, n.pops = 20)
write.csv(results$partition.df, "input/bovis/rhierbaps_clusters_bovis.csv", row.names = FALSE)


###########################
#### BEAST TREE FIGURE ####
###########################

beast_tree <- read.beast("input/bovis/figtree_clusterA_bovis.nexus")

beast_tree <- read.beast("input/bovis/figtree_clusterB_bovis.nexus")

tip_beast <- data.frame(read.csv("input/bovis/metadata_bovis.csv"))

# Merge data with tree
annotated_beast_tree <- full_join(beast_tree, tip_beast, by = c("label" = "label"))

### Check basic tree
ggtree(annotated_beast_tree, mrsd = "2022-11-02", aes(color = Clades)) +
  theme_tree2()


### Base tree
p <- ggtree(annotated_beast_tree, mrsd = "2022-11-02", aes(color = Clades)) +
  ylim(-5, 440) +
  theme_tree2()

### Vertical lines to demarcate sampling periods
p <- p + geom_vline(xintercept = 2004, linetype = "dotted", color = "darkgrey")
p <- p + geom_vline(xintercept = 2015, linetype = "dotted", color = "darkgrey")
p <- p + geom_vline(xintercept = 2022, linetype = "dotted", color = "darkgrey")

p

### Convert heatmap columns to character
tip_beast$IncI     <- as.character(tip_beast$IncI)
tip_beast$IncHI2A  <- as.character(tip_beast$IncHI2A)
tip_beast$IncN     <- as.character(tip_beast$IncN)
tip_beast$rep_1760 <- as.character(tip_beast$rep_1760)


### Build heatmap data frames
heat_IncI     <- data.frame(IncI     = tip_beast$IncI,     row.names = tip_beast$label)
heat_IncHI2A  <- data.frame(IncHI2A  = tip_beast$IncHI2A,  row.names = tip_beast$label)
heat_IncN     <- data.frame(IncN     = tip_beast$IncN,     row.names = tip_beast$label)
heat_rep_1760 <- data.frame(rep_1760 = tip_beast$rep_1760, row.names = tip_beast$label)
heat_Source <- data.frame(Source = tip_beast$Source, row.names = tip_beast$label)

### Numeric column for separator lines
heat_sep <- data.frame(sep = rep(1, nrow(tip_beast)), row.names = tip_beast$label)


### Build composite figure
p2 <- gheatmap(p, heat_IncI, offset = 1.5, width = .02, colnames = TRUE, color = NA,
               font.size = 3, colnames_position = "top", colnames_angle = 45,
               colnames_offset_y = 1, hjust = 0) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                    name = "IncI",
                    labels = c("0" = "Absent", "1" = "Present"),
                    na.value = "grey90") +
  new_scale_fill()
p2

p2 <- gheatmap(p2, heat_IncHI2A, offset = 3.5, width = .02, colnames = TRUE, color = NA,
               font.size = 3, colnames_position = "top", colnames_angle = 45,
               colnames_offset_y = 1, hjust = 0) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                    name = "IncHI2A",
                    labels = c("0" = "Absent", "1" = "Present"),
                    na.value = "grey90") +
  new_scale_fill()


p2 <- gheatmap(p2, heat_IncN, offset = 5.5, width = .02, colnames = TRUE, color = NA,
               font.size = 3, colnames_position = "top", colnames_angle = 45,
               colnames_offset_y = 1, hjust = 0) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                    name = "IncN",
                    labels = c("0" = "Absent", "1" = "Present"),
                    na.value = "grey90") +
  new_scale_fill()
p2

p2 <- gheatmap(p2, heat_rep_1760, offset = 7.5, width = .02, colnames = TRUE, color = NA,
               font.size = 3, colnames_position = "top", colnames_angle = 45,
               colnames_offset_y = 1, hjust = 0) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"),
                    name = "rep_1760",
                    labels = c("0" = "Absent", "1" = "Present"),
                    na.value = "grey90") +
  new_scale_fill()

p2
### Separator lines between columns
p2 <- gheatmap(p2, heat_sep, offset = 2.5, width = 0.001, colnames = FALSE, color = NA, font.size = 0) +
  scale_fill_gradient(low = "grey80", high = "grey80", guide = "none") +
  coord_cartesian(ylim = c(-5, 400)) +
  new_scale_fill()

p2 <- gheatmap(p2, heat_sep, offset = 4.5, width = 0.001, colnames = FALSE, color = NA, font.size = 0) +
  scale_fill_gradient(low = "grey80", high = "grey80", guide = "none") +
  coord_cartesian(ylim = c(-5, 400)) +
  new_scale_fill()

p2 <- gheatmap(p2, heat_sep, offset = 6.5, width = 0.001, colnames = FALSE, color = NA, font.size = 0) +
  scale_fill_gradient(low = "grey80", high = "grey80", guide = "none") +
  coord_cartesian(ylim = c(-5, 400)) +
  new_scale_fill()

p2 <- gheatmap(p2, heat_sep, offset = 8.5, width = 0.001, colnames = FALSE, color = NA, font.size = 0) +
  scale_fill_gradient(low = "grey80", high = "grey80", guide = "none") +
  coord_cartesian(ylim = c(-5, 400)) +
  new_scale_fill()

p2
### add Source heatmap
cols <- c(Feline='#ff7d00', Bovine='#8f00ff', Ovine='#ffdd00', Equine='#01befe', Human='#adff02', Canine='#ff006d', Avian= '#008585',Porcine='#062456', UI_animal='#debcff' )
p2 <- gheatmap(p2, heat_Source, offset=10.5, width=0.02, colnames=T, color=NA, font.size=3, 
               colnames_position="top", colnames_angle=45, colnames_offset_y=1, hjust=0) +
  scale_fill_manual(values=cols, name="Source")
p2

### Also add a legend for the tree branches (Clades colour)
p3 <- p2 + 
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))
p3

ggsave("output/time_tree_bovis_B.png", p3, width=25, height=20, units="cm")

#==============================================================================#
#==========================  PART 02: ZOOMING  ================================#
#==============================================================================#

#### zoom for inset plot

beast_tree <- read.beast("input/bovis/figtree_clusterA_bovis.nexus")

beast_tree <- read.beast("input/bovis/figtree_clusterB_bovis.nexus")
### add data
tip_beast <- data.frame(read.csv("input/bovis/metadata_bovis.csv"))
# Merge data with tree
annotated_beast_tree <- full_join(beast_tree, tip_beast, by = c("label" = "label"))

p <- ggtree(annotated_beast_tree, mrsd = "2022-11-02", aes(color = Clades)) +
  theme_tree2()


p <- p + xlim(1965, 2023)

p <- p + theme(legend.position = "none")

### Demarcate bounds of sampling period
p <- p + geom_vline(xintercept = 2004, linetype = "dotted", color = "darkgrey")
p <- p + geom_vline(xintercept = 2015, linetype = "dotted", color = "darkgrey")
p <- p + geom_vline(xintercept = 2022, linetype = "dotted", color = "darkgrey")

### Read in time series of case data per clade
ts <- read.csv("input/bovis/bovis_time_series.csv")
ts <- as.data.frame(ts)
ts$A <- as.numeric(ts$A)
ts$B <- as.numeric(ts$B)
ts$C <- as.numeric(ts$C)
ts$D <- as.numeric(ts$D)
ts$E <- as.numeric(ts$E)
ts$F <- as.numeric(ts$F)
ts$G <- as.numeric(ts$G)
ts$H <- as.numeric(ts$H)
ts$Date <- as.Date(ts$Date, "%d/%m/%Y")
ts$Date2 <- decimal_date(ts$Date)

# --- Clade colours (update to match your tree) ---
#clade_colours <- c(
#  A = "red",  # red
#  B = "blue",  # blue
#  C = "#99CC66"  # green
#)

clade_colours <- c(
  A = "red",  # red
  B = "gold",  # gold
  C = "#99CC66",  # green
  D = "#4AA155",  # darker green
  E = "#66CCCC",  # orange
  F = "blue",  # brown
  G = "purple",  # pink
  H = "magenta"   # grey
)


# Helper to build each clade inset

make_clade_inset <- function(data, clade, colour, ymax = 2) {
  ggplot(data, aes(x = Date2, y = .data[[clade]], group = 1)) +
    geom_line(colour = colour, linewidth = 0.8) +
    ylab(paste("Clade", clade)) +
    ylim(0, ymax) +
    xlim(1965, 2023) +
    ggfun::theme_noxaxis(axis.line = element_line(linetype = 0)) +
    labs(x = NULL)
}

inset_A <- make_clade_inset(ts, "A", clade_colours["A"])
inset_B <- make_clade_inset(ts, "B", clade_colours["B"])
inset_C <- make_clade_inset(ts, "C", clade_colours["C"])
inset_D <- make_clade_inset(ts, "D", clade_colours["D"])
inset_E <- make_clade_inset(ts, "E", clade_colours["E"])
inset_F <- make_clade_inset(ts, "F", clade_colours["F"])
inset_G <- make_clade_inset(ts, "G", clade_colours["G"])
inset_H <- make_clade_inset(ts, "H", clade_colours["H"])

### Read in effective population size data
tsNe <- read.csv("input/bovis/bovis_popsize_skyline.csv")
tsNe$Ne  <- as.numeric(tsNe$median)
tsNe$LCI <- as.numeric(tsNe$lower)
tsNe$UCI <- as.numeric(tsNe$upper)

inset_Ne <- ggplot(tsNe, aes(x = time, y = Ne)) +
  geom_line(colour = "black", linewidth = 0.8) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "darkgrey", alpha = 0.3) +
  ylab("Ne \n 95% CrI") +
  xlim(1965, 2023) +
  labs(x = NULL)

#==============================================================================#
#=======================  PART 03: ADD TIME SERIES  ===========================#
#==============================================================================#

### Add time series plots to tree
# Height per clade reduced (0.10) to fit 8 clades + Ne without making the figure too tall
p2 <- p %>%
  insert_bottom(inset_A, height = 0.10) %>%
  insert_bottom(inset_B, height = 0.10) %>%
  insert_bottom(inset_C, height = 0.10) %>%
  insert_bottom(inset_D, height = 0.10) %>%
  insert_bottom(inset_E, height = 0.10) %>%
  insert_bottom(inset_F, height = 0.10) %>%
  insert_bottom(inset_G, height = 0.10) %>%
  insert_bottom(inset_H, height = 0.10) %>%
  insert_bottom(inset_Ne, height = 0.25)

p2

# Increased height to accommodate extra panels
ggsave("output/bovis_time_tree_zoom_B.png", p2, width = 25, height = 45, units = "cm")

#==============================================================================#
#=======================  PART 04: ADD OVERLAY GRAPH  =========================#
#==============================================================================#

# Composite image
base_image_path    <- "output/time_tree_bovis.png"
overlay_image_path <- "output/bovis_time_tree_zoom.png"
output_image_path  <- "output/bovis_time_tree_zoom_overlay.png"

base_image_path    <- "output/time_tree_bovis_B.png"
overlay_image_path <- "output/bovis_time_tree_zoom_B.png"
output_image_path  <- "output/bovis_time_tree_zoom_overlay_B.png"

base_image    <- image_read(base_image_path)
overlay_image <- image_read(overlay_image_path)
overlay_image_resized <- image_scale(overlay_image, "1700x1700")
combined_image <- image_composite(base_image, overlay_image_resized, offset = "+50+50")
image_write(combined_image, path = output_image_path, format = "png")
