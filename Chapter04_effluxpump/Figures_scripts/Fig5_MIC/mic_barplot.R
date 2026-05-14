# ============================================================
# MIC Bar Plot - grouped by isolate, 5 bars per isolate
# Aline Rodrigues - Chapter 4/5 biocide & AMR phenotypic data
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 04 Campy cmeABC efflux pump/WRITING/Figures/FIG_MIC/")


# ---- Data entry ----
# Biocide values are in mg/L (converted from %)
# x = not tested (NA)
# S* in result = phenotypically susceptible despite disc diffusion R

data <- tribble(
  ~ID,          ~Species,   ~TET,   ~CIP,   ~ERY,   ~Oxonia,  ~Antigerm,
  "SC0012",     "C. coli",  NA,      NA,   2,      217.5,    9.765,
  "SC0026",     "C. jejuni",32,     4,      NA,     217.5,    9.765,
  "SC0069",     "C. jejuni",64,     8,      NA,     217.5,    19.53,
  "SC0143",     "C. jejuni",NA,   4,      NA,     217.5,    4.88,
  "SC0164",     "C. jejuni",64,     NA,   NA,     217.5,    9.765,
  "SC0169",     "C. jejuni",64,     NA,   NA,     217.5,    9.765,
  "SC0230",     "C. jejuni",32,     4,      NA,     217.5,    9.765,
  "SC0235",     "C. jejuni",32,     4,      NA,     217.5,    9.765,
  "SC0261",     "C. jejuni",32,     4,      NA,     217.5,    4.88,
  "SC0431",     "C. jejuni",32,     4,      NA,     217.5,    9.765,
  "SC0682",     "C. coli",  32,     4,      32,     217.5,    9.765,
  "SC0745",     "C. coli",  32,     8,      NA,     217.5,    19.53,
  "SC0895",     "C. jejuni",NA,   2,      NA,     217.5,    9.765,
  "SC1016",     "C. jejuni",64,     2,      NA,     217.5,    9.765,
  "SC1034",     "C. jejuni",32,     8,      NA,     217.5,    9.765,
  "SC1157",     "C. jejuni",NA,   4,      NA,     217.5,    4.88,
  "SC1357",     "C. jejuni",32,     8,      NA,     217.5,    4.88,
  "SC1513",     "C. jejuni",64,     2,      NA,     217.5,    4.88,
  "SC1704",     "C. jejuni",32,     NA,     NA,     NA,       NA,
  "SC1711",     "C. jejuni",16,     NA,   NA,     217.5,    9.765,
  "NCTC11168",  "C. jejuni",NA,   NA,   NA,     217.5,    4.88,
)

# ---- Reshape to long format ----
data_long <- data %>%
  pivot_longer(
    cols = c(TET, CIP, ERY, Oxonia, Antigerm),
    names_to  = "Drug",
    values_to = "MIC"
  ) %>%
  mutate(
    Drug = factor(Drug,
                  levels = c("TET", "CIP", "ERY", "Oxonia", "Antigerm"),
                  labels = c("Tetracycline",
                             "Ciprofloxacin",
                             "Erythromycin",
                             "Oxonia Active",
                             "Antigerm")),
    ID = factor(ID, levels = unique(data$ID))
  )

# ---- Colour palette (colourblind-friendly) ----
drug_colours <- c(
  "Tetracycline"  = "#E76F51",
  "Ciprofloxacin" = "#2A9D8F",
  "Erythromycin"  = "#457B9D",
  "Oxonia Active" = "#F4A261",
  "Antigerm"      = "#8ECAE6"
)

# ---- Plot ----
p <- ggplot(data_long, aes(x = ID, y = MIC, fill = Drug)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge2(width = 0.9, preserve = "single", padding = 0.05),
    na.rm    = TRUE,
    colour   = "white",
    linewidth = 0.25
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 217.5),
    labels = c("0.01", "0.03", "0.06", "0.12", "0.25", "0.5",
               "1", "2", "4", "8", "16", "32", "64", "128", "217.5"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(values = drug_colours, name = "Drug / Biocide") +
  labs(
    title    = "Minimum Inhibitory Concentrations (MIC) by isolate",
    subtitle = "Bars absent = not tested (NT).",
    x        = "Isolate ID",
    y        = "MIC (mg/L)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    plot.caption     = element_text(size = 8, colour = "grey50", hjust = 0)
  )

print(p)

# ---- Save ----
ggsave(
  "mic_barplot.pdf",
  plot   = p,
  width  = 14,
  height = 6,
  device = cairo_pdf
)


message("Saved: mic_barplot.pdf and mic_barplot.png")
