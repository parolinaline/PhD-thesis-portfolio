
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/GitHub/PhD-thesis-portfolio/Chapter03_Salmonella/Appendix/Apx11_SNPcount_Give/")
getwd()

#Comparing SNPs from Give ST654 from same and different sources

library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Step 1: Read and melt the SNP distance matrix
data <- read.csv("SNP_count_Give_clair.csv")
long <- melt(data)
write.csv(long, "longmatrix_Give_all.csv")

# Step 2: Edit the longmatrix deleting the first column that just had the number of lines 
# and add the ST for both columns with isolates
# Now load the edited file here to filter only the ST654:
data_long <- read.csv("longmatrix_Give_all.csv")
# Filter to keep only ST654 vs ST654 comparisons
data_long <- data_long[data_long$ST1 == 654 & data_long$ST2 == 654, ]
write.csv(data_long, "longmatrix_Give_ST654.csv")

#==============================================================================#
#=====================  IF FILES ARE READY START HERE =========================#
#==============================================================================#
#Comparing isolates from same source and different source

# Step 3: Process to add to "longmatrix_Give_ST654.csv" the period
# For isolates where source1 == source2: same_source = 0
# For isolates where source1 != period2: same_source = 1
data2 <- read.csv("long_Give_ST654_ready.csv")
str(data2)

# Subset by source comparison
long_same <- subset(data2, same_source == 0)
long_diff <- subset(data2, same_source == 1)

# Inspect the distributions
summary(long_same$value)
summary(long_diff$value)

# Compare groups statistically
wilcox.test(long_same$value, long_diff$value)
# p value needs to be lower than 0.05 to be significant
#A p-value < 0.05 means you reject the null hypothesis (that there's no difference between the groups)

# Create a label column for plotting
data2$source_comparison <- ifelse(data2$same_source == 0, "Same source", "Different source")


# ---- Plot 1: Histogram with facets ----
p_hist <- ggplot(data2, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black") +
  facet_wrap(~ source_comparison, ncol = 1, scales = "free_y") +
  labs(x = "Pairwise SNP distance",
       y = "Frequency",
       title = "Pairwise SNP distance Give ST654") +
  theme_bw()
p_hist


# ---- Plot 2: Boxplot with significance annotation ----
p_box <- ggplot(data2, aes(x = source_comparison, y = value, fill = source_comparison)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("Same source", "Different source"))) +
  labs(x = "", 
       y = "Pairwise SNP distance",
       title = "SNP distance by source comparison") +
  scale_fill_manual(values = c("lightgrey", "skyblue")) +
  theme_bw() +
  theme(legend.position = "none")
p_box


# ---- COMBINE PLOTS  ----
p_combined <- p_hist | p_box
p_combined


# ---- Save both plots to PDF ----
pdf("snpdist_Give_source.pdf", width = 15, height = 8, onefile = TRUE)
print(p_combined)
dev.off()
