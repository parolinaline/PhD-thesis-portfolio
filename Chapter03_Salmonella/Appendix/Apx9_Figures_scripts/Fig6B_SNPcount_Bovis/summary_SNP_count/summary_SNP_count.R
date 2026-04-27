
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/6D_SNPcount_Bovis/summary_SNP_count/")

getwd()

################################################################################
###################### ALL ST377 2003-2022 #####################################
################################################################################


###########################
### INPUT FILES ###########
###########################
# Pairwise SNP matrix (square matrix of SNP counts)
snps <- read.csv("Bovis_ST377_SNPs.csv", row.names = 1, check.names = FALSE)

# Metadata file (must have columns: isolate, source)
meta <- read.csv("metadata_ST377.csv")

###########################
### PREPARE DATA ##########
###########################
# Melt the SNP matrix into long format
snps_long <- melt(as.matrix(snps), varnames = c("isolate1", "isolate2"), value.name = "SNPs")

# Remove self-comparisons (distance = 0)
snps_long <- snps_long %>% filter(isolate1 != isolate2)

# Join source metadata for both isolates
snps_long <- snps_long %>%
  left_join(meta, by = c("isolate1" = "ID")) %>%
  rename(source1 = Source) %>%
  left_join(meta, by = c("isolate2" = "ID")) %>%
  rename(source2 = Source)

###########################
### CLASSIFY COMPARISONS ##
###########################
# Define same vs different source comparisons
snps_long <- snps_long %>%
  mutate(comparison_type = ifelse(source1 == source2, "Same source", "Different sources"))

###########################
### SUMMARIZE STATS #######
###########################
# Overall SNP summary
overall_summary <- snps_long %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs)
  )

# Summary by source (median and max SNPs among isolates from the same source)
same_source_summary <- snps_long %>%
  filter(comparison_type == "Same source") %>%
  group_by(source1) %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs),
    .groups = "drop"
  )

# Summary for comparisons across sources
diff_source_summary <- snps_long %>%
  filter(comparison_type == "Different sources") %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs)
  )

###########################
### PRINT RESULTS #########
###########################
cat("\n===== Overall SNP Summary =====\n")
print(overall_summary)

cat("\n===== Same-source SNP Summary =====\n")
print(same_source_summary)

cat("\n===== Different-source SNP Summary =====\n")
print(diff_source_summary)


################################################################################
################### ALL ST377 ISOLATED AFTER 2015 ##############################
################################################################################


###########################
### LOAD PACKAGES #########
###########################
library(tidyverse)
library(reshape2)

###########################
### INPUT FILES ###########
###########################
# Pairwise SNP matrix (square matrix of SNP counts)
snps <- read.csv("Bovis_ST377_2015to2022.csv", row.names = 1, check.names = FALSE)

# Metadata file (must have columns: isolate, source)
meta <- read.csv("metadata_2015to2022.csv")

###########################
### PREPARE DATA ##########
###########################
# Melt the SNP matrix into long format
snps_long <- melt(as.matrix(snps), varnames = c("isolate1", "isolate2"), value.name = "SNPs")

# Remove self-comparisons (distance = 0)
snps_long <- snps_long %>% filter(isolate1 != isolate2)

# Join source metadata for both isolates
snps_long <- snps_long %>%
  left_join(meta, by = c("isolate1" = "ID")) %>%
  rename(source1 = Source) %>%
  left_join(meta, by = c("isolate2" = "ID")) %>%
  rename(source2 = Source)

###########################
### CLASSIFY COMPARISONS ##
###########################
# Define same vs different source comparisons
snps_long <- snps_long %>%
  mutate(comparison_type = ifelse(source1 == source2, "Same source", "Different sources"))

###########################
### SUMMARIZE STATS #######
###########################
# Overall SNP summary
overall_summary <- snps_long %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs)
  )

# Summary by source (median and max SNPs among isolates from the same source)
same_source_summary <- snps_long %>%
  filter(comparison_type == "Same source") %>%
  group_by(source1) %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs),
    .groups = "drop"
  )

# Summary for comparisons across sources
diff_source_summary <- snps_long %>%
  filter(comparison_type == "Different sources") %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs)
  )

###########################
### PRINT RESULTS #########
###########################
cat("\n===== Overall SNP Summary =====\n")
print(overall_summary)

cat("\n===== Same-source SNP Summary =====\n")
print(same_source_summary)

cat("\n===== Different-source SNP Summary =====\n")
print(diff_source_summary)


################################################################################
############ ALL ST377 ISOLATED AFTER 2015 EXCLUDING SBV313 ####################
################################################################################

###########################
### INPUT FILES ###########
###########################
# Pairwise SNP matrix (square matrix of SNP counts)
snps <- read.csv("Bovis_ST377_recent_no313.csv", row.names = 1, check.names = FALSE)

# Metadata file (must have columns: isolate, source)
meta <- read.csv("metadata_recent_no313.csv")

###########################
### PREPARE DATA ##########
###########################
# Melt the SNP matrix into long format
snps_long <- melt(as.matrix(snps), varnames = c("isolate1", "isolate2"), value.name = "SNPs")

# Remove self-comparisons (distance = 0)
snps_long <- snps_long %>% filter(isolate1 != isolate2)

# Join source metadata for both isolates
snps_long <- snps_long %>%
  left_join(meta, by = c("isolate1" = "ID")) %>%
  rename(source1 = Source) %>%
  left_join(meta, by = c("isolate2" = "ID")) %>%
  rename(source2 = Source)

###########################
### CLASSIFY COMPARISONS ##
###########################
# Define same vs different source comparisons
snps_long <- snps_long %>%
  mutate(comparison_type = ifelse(source1 == source2, "Same source", "Different sources"))

###########################
### SUMMARIZE STATS #######
###########################
# Overall SNP summary
overall_summary <- snps_long %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs)
  )

# Summary by source (median and max SNPs among isolates from the same source)
same_source_summary <- snps_long %>%
  filter(comparison_type == "Same source") %>%
  group_by(source1) %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs),
    .groups = "drop"
  )

# Summary for comparisons across sources
diff_source_summary <- snps_long %>%
  filter(comparison_type == "Different sources") %>%
  summarise(
    median_SNPs = median(SNPs),
    max_SNPs = max(SNPs)
  )

###########################
### PRINT RESULTS #########
###########################
cat("\n===== Overall SNP Summary =====\n")
print(overall_summary)

cat("\n===== Same-source SNP Summary =====\n")
print(same_source_summary)

cat("\n===== Different-source SNP Summary =====\n")
print(diff_source_summary)
