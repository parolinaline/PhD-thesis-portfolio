## Description
Figure 2  - Distribution of plasmid  marker combinations per plasmid contig and plasmid contig lengths across Campylobacter STs. Heatmap showing the distribution of plasmid marker combinations (columns) across STs (rows) for C. jejuni (A) and C. coli (B). The colour gradient represents the number of isolates within each ST that carry a given plasmid marker combination. To the right of the heatmap, a horizontal boxplot displays the distribution of plasmid contig lengths for each ST in C. jejuni (A) and C. coli (B). Individual dots represent single plasmid contigs, providing a visual estimate of the number of plasmids per ST. Figures were generated using ggplot2 and arranged with aplot in R.

### *Campylobacter jejuni* isolates

#### Input files
[plasmid markers for *jejuni*](./jejuni_markers.csv)
[boxplot](./boxplot_jejuni.csv)
[output](./fig02_j.pdf)


```r

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/Writing_Chapter05/Thesis_Figures_scripts_data/Fig02_ilm_plasmid")

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(aplot)

# Read your Excel or CSV file
data <- read.csv("jejuni_markers.csv")  # replace with readxl::read_excel() if needed

# Count number of isolates per ST that have each marker
count_data <- data %>%
  distinct(ID, GENES, ST) %>%  # prevent duplicate counts per isolate
  group_by(ST, GENES) %>%
  summarise(N = n_distinct(ID), .groups = 'drop')

# Pivot to wide format
wide_data <- count_data %>%
  pivot_wider(names_from = GENES, values_from = N, values_fill = 0)

# Set ST as rownames
heatmap_matrix <- as.data.frame(wide_data)
rownames(heatmap_matrix) <- heatmap_matrix$ST
heatmap_matrix$ST <- NULL

#Melt Data for ggplot
heatmap_long <- melt(as.matrix(heatmap_matrix))
colnames(heatmap_long) <- c("ST", "Plasmid_marker", "Count")

# Vector these STs to be highlightes (they were chosen to be Nanopored)
#highlight_STs <- c("ST422","ST51","ST704","ST14186","ST2357","ST520","ST53","ST42","ST3610","ST1030","ST22","ST45","ST48","ST677","ST50",
#                   "ST61","ST12762","ST2392","ST5647","ST137","ST535")

#Create a colour vector for all STs
#row_colors <- ifelse(levels(factor(heatmap_long$ST)) %in% highlight_STs, "red", "black")
#names(row_colors) <- levels(factor(heatmap_long$ST))


#Plot the Heatmap
heatmap_jejuni<- ggplot(heatmap_long, aes(x = Plasmid_marker, y = ST, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "#191970", name = "Number of isolates") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
#    ,
#    axis.text.y = element_text(color = row_colors)
  ) +
  ylab("Sequence Type (ST)") +
  xlab(NULL) +
  ggtitle(expression("A) Plasmid marker combinations per plasmid contig in "*italic("C. jejuni")))

heatmap_jejuni

##########################BOXPLOT##############################################
########################PLASMID SIZE############################################
################################################################################
#Load data
boxplot_jejuni <- read.csv("boxplot_jejuni.csv")


boxplot_figure_jejuni <- ggplot(boxplot_jejuni, aes(y = ST, x = Length)) +
  geom_boxplot(fill = "lightblue", color = "darkgrey") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal() +
  labs(x = "Distribution of plasmid contig lengths per ST") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

boxplot_figure_jejuni


#### PUT TOGETHER BOXPLOT AND HEATMAP

final_plot_jejuni <- heatmap_jejuni %>%
  insert_right(boxplot_figure_jejuni, width = 0.8)


final_plot_jejuni
################################################################################
################################################################################
# Save the figure as SVG
#install.packages("svglite")
library(svglite)
ggsave("ILMN_jejuni.svg", plot = final_plot_jejuni, width = 18, height = 12, units = "in", dpi = 600)

```

### *Campylobacter coli* isolates

#### Input files
[plasmid markers for *coli*](./coli_markers.csv)
[boxplot](./boxplot_coli.csv)
[output](./fig02_c.pdf)

```r

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/Writing_Chapter05/Thesis_Figures_scripts_data/Fig02_ilm_plasmid")

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(aplot)

# Read your Excel or CSV file
data_coli <- read.csv("coli_markers.csv")  # replace with readxl::read_excel() if needed

# Count number of isolates per ST that have each marker
count_data_coli <- data_coli %>%
  distinct(ID, GENES, ST) %>%  # prevent duplicate counts per isolate
  group_by(ST, GENES) %>%
  summarise(N = n_distinct(ID), .groups = 'drop')

# Pivot to wide format
wide_data_coli <- count_data_coli %>%
  pivot_wider(names_from = GENES, values_from = N, values_fill = 0)

# Set ST as rownames
heatmap_matrix_coli <- as.data.frame(wide_data_coli)
rownames(heatmap_matrix_coli) <- heatmap_matrix_coli$ST
heatmap_matrix_coli$ST <- NULL

#Melt Data for ggplot
heatmap_long_coli <- melt(as.matrix(heatmap_matrix_coli))
colnames(heatmap_long_coli) <- c("ST", "Plasmid_marker", "Count")

# Vector these STs to be highlightes (they were chosen to be Nanopored)
#highlight_STs_coli <- c("ST2256","ST8926","ST900","ST4009","ST14595","ST830","ST3230","ST829","ST1115","ST3232","ST3222","ST1096")

#Create a colour vector for all STs
#row_colors_coli <- ifelse(levels(factor(heatmap_long_coli$ST)) %in% highlight_STs_coli, "red", "black")
#names(row_colors_coli) <- levels(factor(heatmap_long_coli$ST))


#Plot the Heatmap
heatmap_coli <- ggplot(heatmap_long_coli, aes(x = Plasmid_marker, y = ST, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightyellow", high = "#191970", name = "Number of isolates") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
#    ,
#    axis.text.y = element_text(color = row_colors_coli)
  ) +
  ylab("Sequence Type (ST)") +
  xlab(NULL) +
  ggtitle(expression("B) Plasmid marker combinations per contig in "*italic("C. coli")))

heatmap_coli

##########################BOXPLOT##############################################
########################PLASMID SIZE############################################
################################################################################
#Load data
boxplot_coli <- read.csv("boxplot_coli.csv")


boxplot_figure_coli <- ggplot(boxplot_coli, aes(y = ST, x = Length)) +
  geom_boxplot(fill = "lightblue", color = "darkgrey") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal() +
  labs(x = "Distribution of plasmid contig lengths per ST") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())

boxplot_figure_coli


#### PUT TOGETHER BOXPLOT AND HEATMAP

final_plot_coli <- heatmap_coli %>%
  insert_right(boxplot_figure_coli, width = 0.8)


final_plot_coli
################################################################################
################################################################################
# Save the figure as SVG
#install.packages("svglite")
library(svglite)
ggsave("ilmn_marker_coli2.svg", plot = final_plot_coli, width = 15, height = 10, units = "in", dpi = 600)

```
