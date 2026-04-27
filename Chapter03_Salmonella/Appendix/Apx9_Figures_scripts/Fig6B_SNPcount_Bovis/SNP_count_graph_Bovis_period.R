#Comparing SNPs from Bovis ST377 but isolates in two different periods:
# period 01: before 2015
# period 02: from 2015 onwards
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/6D_SNPcount_Bovis/")

library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Step 1: Read and melt the SNP distance matrix
data <- read.csv("SNPcount_Bovis_clair.csv")
long <- melt(data)
write.csv(long, "longmatrix_Bovis_all.csv")

# Step 2: Edit the longmatrix deleting the first column that just had the number of lines 
# and add the ST for both columns with isolates
# Now load the edited file here to filter only the ST377:
data_long <- read.csv("longmatrix_Bovis_all.csv")
# Filter to keep only ST377 vs ST377 comparisons
data_long <- data_long[data_long$ST1 == 377 & data_long$ST2 == 377, ]
write.csv(data_long, "longmatrix_Bovis_ST377.csv")

#==============================================================================#
#=====================  IF FILES ARE READY START HERE =========================#
#==============================================================================#

# Step 3: Process to add to "longmatrix_Bovis_ST377.csv" the year of isolation from each 
# of the ID columns and name them year1 and year2
# 2004 until 2014 are named period01 and 2015 to 2022 are named period02
# For isolates where year1 == year2: same_period = 0
# For isolates where year1 != year2: same_period = 1
data2 <- read.csv("long_Bovis_ST377_ready.csv")
str(data2)

# Subset by period comparison
long_same <- subset(data2, same_period == 0)
long_diff <- subset(data2, same_period == 1)

# Inspect the distributions
summary(long_same$value)
summary(long_diff$value)

# Compare groups statistically
wilcox.test(long_same$value, long_diff$value)
# p value needs to be lower than 0.05 to be significant
#A p-value < 0.05 means you reject the null hypothesis (that there's no difference between the groups)

# Create a label column for plotting
data2$period_comparison <- ifelse(data2$same_period == 0, "Same period", "Different period")

# ---- Plot 1: Histogram with facets ----
p_hist <- ggplot(data2, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black") +
  facet_wrap(~ period_comparison, ncol = 1, scales = "free_y") +
  labs(x = "Pairwise SNP distance",
       y = "Frequency",
       title = "Pairwise SNP distance Bovismorbificans ST377") +
  theme_bw()
p_hist


# ---- Plot 2: Boxplot with significance annotation ----
p_box <- ggplot(data2, aes(x = period_comparison, y = value, fill = period_comparison)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("Same period", "Different period"))) +
  labs(x = "", 
       y = "Pairwise SNP distance",
       title = "SNP distance by period comparison") +
  scale_fill_manual(values = c("lightgrey", "skyblue")) +
  theme_bw() +
  theme(legend.position = "none")
p_box


# ---- COMBINE PLOTS  ----
p_combined <- p_hist | p_box
p_combined


# ---- Save both plots to PDF ----
pdf("snpdist_Bovis_periodnew.pdf", width = 15, height = 8, onefile = TRUE)
print(p_combined)
dev.off()
