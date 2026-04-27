
getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/Fig_TET_resistance/")

# 📦 Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# 📂 Read your Excel file
df <- read_excel("interpreted_data_upsetplot.xlsx", sheet = "Bovis")

# 🧮 Summarise per Year and ST
summary_df <- df %>%
  group_by(Year, ST) %>%
  summarise(
    Total = n(),
    Resistant = sum(TET == "R", na.rm = TRUE),
    .groups = "drop"
  )

# 🖼 Panel A: total isolates per year (all STs combined) + labels
totals_df <- summary_df %>%
  group_by(Year) %>%
  summarise(Total = sum(Total))

p1 <- ggplot(totals_df, aes(x = factor(Year), y = Total)) +
  geom_line(group = 1, color = "black") +  # group=1 needed for factor x
  geom_point(color = "black", size = 3) +
  geom_text(aes(label = Total), vjust = -0.8, size = 3.5) +
  scale_x_discrete(drop = FALSE) +  # ensures all levels show
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    axis.line = element_line()
  ) +
  labs(
    title = "A. Total isolates per year",
    y = "Number of isolates",
    x = NULL
  )

p1


# 🖼 Panel B: heatmap (Year × ST) of resistant counts
p2 <- ggplot(summary_df, aes(x = factor(Year), y = factor(ST), fill = Resistant)) +
  geom_tile(color = NA) +  # no borders on the tiles
  geom_text(data = subset(summary_df, Resistant > 0),
            aes(label = Resistant), color = "black", size = 3.5) +
  scale_fill_gradient(low = "white", high = "red", name = "No. isolates") +
#  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "B. Tetracycline-resistant isolates per ST",
    x = "Year",
    y = "Sequence type"
  )
p2

# Combine plots
p1 / p2 + plot_layout(heights = c(0.8, 0.2))


ggsave("heatmap_lineplot_tet_Bovis.pdf", width = 10, height = 7)
