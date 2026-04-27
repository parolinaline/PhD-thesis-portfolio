# Install necessary packages if you haven't already
# Load the packages
library(ggplot2)
library(RColorBrewer)
library(ComplexUpset)
library(dplyr)
library(readxl)


getwd()

setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/5A_UpsetPlot_Plasmid_CoExistence")

#SOURCE OF THE CODE:
#https://krassowski.github.io/complex-upset/articles/Examples_R.html


################################################################################
########### UPSET PLOT - PLASMID INCOMPATIBILITY GROUPS ########################
################################################################################

# Load your binary plasmid data
# The file should have:
# Column 1 = ID
# Column 2 = ST
# Columns 3+ = incompatibility groups (0/1 matrix)


plasmid_data <- read_excel("Matrix_Inc_Bovis.xlsx", sheet = "All")

# Separate metadata and binary matrix
metadata <- plasmid_data %>% select(1:2)  # ID and ST
binary_matrix <- plasmid_data %>% select(-c(1:2))  # all incompatibility group columns

# Combine metadata and binary columns
binary_data <- bind_cols(metadata, binary_matrix)

# Check structure
str(binary_data)

# Extract Inc group names automatically
inc_groups <- colnames(binary_matrix)
names(inc_groups) <- inc_groups

################################################################################
########### BUILD UPSET PLOT ###################################################
################################################################################

upset(binary_data,
      inc_groups,
      name = "Plasmid type combinations",
      set_sizes = (
        upset_set_size() +
          ylab("No. of isolates")
      ),
      base_annotations = list(
        'Intersection size' = intersection_size(
          counts = TRUE,
          mapping = aes(fill = ST),
          text = list(size = 6, color = "black")
        ) +
          geom_bar(stat = "identity") +
          ylab("Number of isolates") +
          labs(fill = "Sequence type (ST)") +
          ggtitle(expression(paste("Co-existence of plasmid types in ", italic("Salmonella"), " Bovismorbificans"))) +
          theme(
            axis.title.y = element_text(size = 18, face = "bold"),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14),
            plot.title = element_text(size = 20, face = "bold")
          )
      ),
      themes = upset_modify_themes(list(
        'intersections_matrix' = theme(text = element_text(size = 12, face = "bold")),
        'overall_sizes' = theme(axis.text.x = element_text(size = 12))
      )),
      matrix = (
        intersection_matrix(
          geom = geom_point(shape = 'square', size = 3),
          segment = geom_segment(linewidth = 1.3)
        )
      ),
      width_ratio = 0.1,
      stripes = upset_stripes(
        geom = geom_segment(size = 8),
        colors = c('#E5E4E2', "#fdfdfd")
      ),
      sort_sets = FALSE
) & scale_fill_manual(values = brewer.pal(8, "Set1"))

# Save output
ggsave("Upset_Plasmid_IncGroups.pdf", width = 18, height = 15)

