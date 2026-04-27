# Install necessary packages if you haven't already
# Load the packages
library(ggplot2)
library(RColorBrewer)
library(ComplexUpset)
library(dplyr)
library(readxl)

install.packages('ComplexUpset')
#remotes::install_version("ggplot2", version = "3.5.1")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("krassowski/complex-upset")


getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/2_Fig_upsetplot_AMR/")


#SOURCE OF THE CODE:
#https://krassowski.github.io/complex-upset/articles/Examples_R.html

################################################################################
###########BOVISMORBIFICANS BY SOURCE###########################################

# Load your data
# Replace 'Sheet1' with the actual sheet name if it's different


#USING CLSI INTERPRETATION
data_bovis <- read_excel("interpreted_data_upsetplot.xlsx", sheet = "Bovis")

# Select metadata columns (1-5) and result columns (6-19)
metadata_bovis <- data_bovis %>% select(1:4)
results_bovis <- data_bovis %>% select(5:18)

# Convert result columns from R, S, I to binary (1 for R, 0 for S and I)
binary_results_Bovis <- results_bovis %>%
  mutate(across(everything(), ~ ifelse(. == "R", 1, 0)))


# Structure of the variables, seeing that there is the binary 0,1 matrix
str(binary_results_Bovis)

# Combine the metadata and binary results
binary_data_Bovis <- bind_cols(metadata_bovis, binary_results_Bovis)

# Specifying names used for the antibiotic labels
antibiotics <- colnames(results_bovis)  # Automatically get column names for antibiotics
names(antibiotics) <- antibiotics

# Constructing the plot with ComplexUpset
PLOT <- upset(binary_data_Bovis,
      antibiotics,
      name = "Resistance combinations",
      set_sizes = ( upset_set_size()+
                      ylab("No. of isolates")),
      themes = upset_modify_themes(
        list(
          'intersections_matrix'= theme(text=element_text(size = 20, face = "bold")),
          'overall_sizes'=theme(axis.text.x = element_text(size = 12)))),
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts =TRUE,
          mapping=aes(fill = Source),
          text = list(size=8, color ="black")
        ) 
        + geom_bar(stat="identity")
        + ylab('Number of isolates')
        + labs(fill = "Source") +
          ggtitle(expression(paste("AMR profile of ", italic("Salmonella"), " Bovismorbificans"))) +
          scale_y_continuous(breaks = seq(0, 340, by = 20)) +
          theme(axis.title.y = element_text(size = 20, face = "bold"),
                axis.text.y = element_text(size = 12),
                legend.title = element_text(size = 18, face = "bold"),
                legend.text = element_text(size = 18),
                plot.title = element_text(size = 22, face = "bold")
          )),
      width_ratio=0.1,
      matrix = (intersection_matrix(geom = geom_point(shape='square', size=4.5),
                                    segment = geom_segment(linewidth = 1.5))),
      sort_sets=FALSE,
      stripes=upset_stripes(geom = geom_segment(size=8),
                            colors=c('#E5E4E2', "#fdfdfd")),
) & scale_fill_manual(values = c( "#ffd3b6","#a786dc","#d2c1ed","#5499c7","#80dbe3","#87c6a2","#0aae24", "#FFD700", "red"))

PLOT


ggsave("Upset_Bovismorbificans2026.pdf", width = 18, height = 15)



################################################################################
#############GIVE BY SOURCE#####################################################

# Load your data
# Replace 'Sheet1' with the actual sheet name if it's different


#USING CLSI INTERPRETATION
data_give <- read_excel("interpreted_data_CLSI.xlsx", sheet = "Give")


# Select metadata columns (1-5) and result columns (6-19)
metadata_give <- data_give %>% select(1:4)
results_give <- data_give %>% select(5:18)

# Convert result columns from R, S, I to binary (1 for R, 0 for S and I)
binary_results_give <- results_give %>%
  mutate(across(everything(), ~ ifelse(. == "R", 1, 0)))


# Structure of the variables, seeing that there is the binary 0,1 matrix
str(binary_results_give)

# Combine the metadata and binary results
binary_data_give <- bind_cols(metadata_give, binary_results_give)

# Specifying names used for the antibiotic labels
antibiotics <- colnames(results_give)  # Automatically get column names for antibiotics
names(antibiotics) <- antibiotics


####UPSET
upset(binary_data_give,
      antibiotics,
      name = "Resistance combinations",
      set_sizes = ( upset_set_size()+
                      ylab("No. of isolates")),
      themes = upset_modify_themes(
        list(
          'intersections_matrix'= theme(text=element_text(size = 20, face = "bold")),
          'overall_sizes'=theme(axis.text.x = element_text(size = 12)))),
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts =TRUE,
          mapping=aes(fill = Source),
          text = list(size=8, color ="black")
        ) 
        + geom_bar(stat="identity") 
        + ylab('Number of isolates')
        + labs(fill = "Source") +
          ggtitle(expression(paste("AMR profile of ", italic("Salmonella"), " Give"))) +
          scale_y_continuous(breaks = seq(0, 80, by = 5)) +
          theme(axis.title.y = element_text(size = 20, face = "bold"),
                axis.text.y = element_text(size = 12),
                legend.title = element_text(size = 18, face = "bold"),
                legend.text = element_text(size = 18),
                plot.title = element_text(size = 22, face = "bold")
          )),
      width_ratio=0.1,
      matrix = (intersection_matrix(geom = geom_point(shape='square', size=4.5),
                                    segment = geom_segment(linewidth = 1.5))),
      sort_sets=FALSE,
      stripes=upset_stripes(geom = geom_segment(size=8),
                            colors=c('#E5E4E2', "#fdfdfd"))
) & scale_fill_manual(values = c("#ffd3b6","#a786dc","#d2c1ed","#5499c7","#80dbe3","#87c6a2","#0aae24"))


ggsave("Upset_Give_CLSI.pdf", width = 15, height = 12)



################################################################################
################# RESISTANCE BY SEROVAR ########################################

# Load your data
# Replace 'Sheet1' with the actual sheet name if it's different
data_both_byserovar <- read_excel("AST_upsetplot_BovisandGive.xlsx")

# Select metadata columns (1-5) and result columns (6-19)
metadata_both_byserovar <- data_both_byserovar %>% select(1:5)
results_both_byserovar <- data_both_byserovar %>% select(6:19)

# Convert result columns from R, S, I to binary (1 for R, 0 for S and I)
binary_results_both_byserovar <- results_both_byserovar %>%
  mutate(across(everything(), ~ ifelse(. == "R", 1, 0)))


# Structure of the variables, seeing that there is the binary 0,1 matrix
str(binary_results_both_byserovar)

# Combine the metadata and binary results
binary_data_both_byserovar <- bind_cols(metadata_both_byserovar, binary_results_both_byserovar)

# Specifying names used for the antibiotic labels
antibiotics <- colnames(results_both_byserovar)  # Automatically get column names for antibiotics
names(antibiotics) <- antibiotics

# Constructing the plot with ComplexUpset
upset(binary_data_both_byserovar,
      antibiotics,
      name = "Resistance combinations for Bovismorbificans and Give", 
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill = Serovar)
        ) 
        + geom_bar(stat="identity") 
        + ylab('Number of isolates')
        + labs(fill = "Serovar") +
          ggtitle(expression(paste("AMR profile of ", italic("Salmonella"), " Bovismorbificans and Give"))) +
          scale_y_continuous(breaks = seq(0, 600, by = 20))),
      width_ratio=0.1,
      stripes='white'
) & scale_fill_manual(values = c("#d2c1ed","#a9dfbf"))


################################################################################
##################### RESISTANCE BY YEAR #######################################

# Load your data
# Replace 'Sheet1' with the actual sheet name if it's different
data_both_byserovar <- read_excel("AST_upsetplot_BovisandGive.xlsx")

# Select metadata columns (1-5) and result columns (6-19)
metadata_both_byserovar <- data_both_byserovar %>% select(1:5)
results_both_byserovar <- data_both_byserovar %>% select(6:19)

# Convert result columns from R, S, I to binary (1 for R, 0 for S and I)
binary_results_both_byserovar <- results_both_byserovar %>%
  mutate(across(everything(), ~ ifelse(. == "R", 1, 0)))


# Structure of the variables, seeing that there is the binary 0,1 matrix
str(binary_results_both_byserovar)

# Combine the metadata and binary results
binary_data_both_byserovar <- bind_cols(metadata_both_byserovar, binary_results_both_byserovar)

# Specifying names used for the antibiotic labels
antibiotics <- colnames(results_both_byserovar)  # Automatically get column names for antibiotics
names(antibiotics) <- antibiotics

# colour pallete

color_palette <- c("#201923", "#ffffff", "#fcff5d", "#7dfc00", "#0ec434", "#228c68", "#8ad8e8", "#235b54", "#29bdab", "#3998f5", 
                   "#37294f", "#7f00ff", "#ff0054", "#ff5400", "#ffbd00", "#6e7f80", "#00c6c1", "#264653", "#2a9d8f", "#e9c46a", 
                   "#f4a261", "#e76f51", "#d4a5a5", "#a786dc", "#ef476f", "#118ab2", "#06d6a0", "#ffd166", "#073b4c", "#8e44ad")

# Constructing the plot with ComplexUpset
upset(binary_data_both_byserovar,
      antibiotics,
      name = "Resistance combinations for Bovismorbificans and Give", 
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill = Year)
        ) 
        + geom_bar(stat="identity") 
        + ylab('Number of isolates')
        + labs(fill = "Year") +
          ggtitle(expression(paste("AMR profile of ", italic("Salmonella"), " Bovismorbificans and Give"))) +
          scale_y_continuous(breaks = seq(0, 600, by = 20))),
      width_ratio=0.1,
      stripes='white'
) & scale_fill_manual(values = color_palette)


################################################################################
##################### INTERMEDIATE RESISTANCE BY SEROVAR########################

# Load your data
# Replace 'Sheet1' with the actual sheet name if it's different
data_int <- read_excel("AST_upsetplot_BovisandGive.xlsx")

# Select metadata columns (1-5) and result columns (6-19)
metadata_int <- data_int %>% select(1:5)
results_int <- data_int %>% select(6:19)

# Convert result columns from R, S, I to binary (1 for R, 0 for S and I)
binary_results_int <- results_int %>%
  mutate(across(everything(), ~ ifelse(. == "I", 1, 0)))


# Structure of the variables, seeing that there is the binary 0,1 matrix
str(binary_results_int)

# Combine the metadata and binary results
binary_data_int <- bind_cols(metadata_int, binary_results_int)

# Specifying names used for the antibiotic labels
antibiotics <- colnames(results_int)  # Automatically get column names for antibiotics
names(antibiotics) <- antibiotics

# Constructing the plot with ComplexUpset
upset(binary_data_int,
      antibiotics,
      name = "Intermediate resistance combinations for Bovismorbificans and Give", 
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill = Serovar)
        ) 
        + geom_bar(stat="identity") 
        + ylab('Number of isolates')
        + labs(fill = "Serovar") +
          ggtitle(expression(paste("Intermediate resistance of ", italic("Salmonella"), " Bovismorbificans and Give"))) +
          scale_y_continuous(breaks = seq(0, 360, by = 20))),
      width_ratio=0.1,
      stripes='white'
) & scale_fill_manual(values = c("#a786dc","#a9dfbf"))
