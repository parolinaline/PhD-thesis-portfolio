# Figure 01

### Description

Decision tree  illustrating the criteria used to select the assembler tool for plasmid assembly. Hybrid assemblies (Unicycler) and long read assemblies (Flye) were compared, and the assembler was chosen based on plasmid detection and the number of circularisable plasmids. Figure generated using DiagrammeR v1.0.11

#### Output
[Output](./decision_tree2.pdf)

```R

# Load and/or install packages

install.packages("DiagrammeR")

library(DiagrammeRsvg)
library(DiagrammeR)
library(rsvg)

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/Writing_Chapter05/Thesis_Figures_scripts_data/Fig01_decisiontree_assembly")

diagram <- grViz("
digraph decision_tree {
  graph [layout = dot, rankdir = TB]

  # --- Nodes ---
  short_read_decision [label = 'Plasmids detected\\nusing short reads only?', shape=diamond, style=filled, fillcolor='#fff9c4']

  sr_yes [label = 'Yes', shape=plaintext]
  sr_no  [label = 'No', shape=plaintext]

  selected64 [label = '64 representative isolates\\nselected for Nanopore sequencing', shape=box, style=filled, fillcolor='#e0f7fa']

  testing_longreads [label = 'Testing assembly tools for long reads:\\nFlye vs Unicycler', shape=box, style=filled, fillcolor='#ffe0b2']

  flye_only_long [label = 'Plasmids detected only \\nin Flye assembly', shape=box, style=filled, fillcolor='#ffccbc']
  uni_only_long  [label = 'Plasmids detected only \\nin Unicycler assembly', shape=box, style=filled, fillcolor='#ffccbc']
  both_long      [label = 'Plasmids detected \\nwith both tools', shape=box, style=filled, fillcolor='#ffccbc']
  none_long      [label = 'No plasmids detected \\nusing long reads', shape=box, style=filled, fillcolor='#ffccbc']

  bandage_count [label = 'Number of complete plasmids\\n(Bandage output)', shape=box, style=filled, fillcolor='#bbdefb']

  # Numbered comparison boxes (connected from bandage_count)
  num1 [label = 'Flye > Unicycler', shape=box, style=filled, fillcolor='#d1c4e9']
  num2 [label = 'Unicycler > Flye', shape=box, style=filled, fillcolor='#d1c4e9']
  num3 [label = 'Equal number of \\ncomplete plasmids', shape=box, style=filled, fillcolor='#d1c4e9']

  # Small follow-up nodes that state chosen tool
  chosen_flye [label = 'Chosen tool:\\nFlye', shape=oval, style=filled, fillcolor='#c8e6c9']
  chosen_unicycler [label = 'Chosen tool:\\nUnicycler', shape=oval, style=filled, fillcolor='#c8e6c9']

  # Excluded node for long-read negative
  excluded_from_analysis [label = 'Excluded from analysis', shape=box, style=filled, fillcolor='#ffcdd2']

  # --- Edges ---
  short_read_decision -> sr_yes
  short_read_decision -> sr_no

  sr_yes -> selected64
  sr_no  -> excluded_from_analysis [label='Not selected for Nanopore', fontsize=10]

  selected64 -> testing_longreads

  testing_longreads -> flye_only_long
  testing_longreads -> uni_only_long
  testing_longreads -> both_long
  testing_longreads -> none_long

  # Direct arrows from single-tool detections to chosen-tool outcomes
  flye_only_long -> chosen_flye
  uni_only_long  -> chosen_unicycler

  # If no plasmids found in long reads -> excluded
  none_long -> excluded_from_analysis

  # From 'both' to Bandage count, then to the numbered outcomes
  both_long -> bandage_count

  bandage_count -> num1
  bandage_count -> num2
  bandage_count -> num3

  # Numbered outcomes link to chosen-tool boxes
  num1 -> chosen_flye
  num2 -> chosen_unicycler
  num3 -> chosen_unicycler

  # --- layout control ---
  { rank = same; flye_only_long; uni_only_long; both_long; none_long }    # keep the 4 long-read outcomes aligned
  { rank = same; num1; num2; num3 }                                      # keep the 3 comparisons aligned
  { rank = same; chosen_flye; chosen_unicycler }                         # final choices aligned horizontally

  # Force 'testing_longreads' directly under selected64
  selected64 -> testing_longreads [minlen=1]

  edge [minlen=1]
}
")

# Convert DiagrammeR output to SVG
svg_code <- DiagrammeRsvg::export_svg(diagram)

# Save as a high-resolution PDF (vector — infinite resolution)
rsvg::rsvg_pdf(
  charToRaw(svg_code),
  file = "C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/decision_tree2.pdf",
  width = 20 * 600,   # 7 inches × 72 dpi
  height = 10 * 600  # 10 inches × 72 dpi (adjust as needed)
)

```
