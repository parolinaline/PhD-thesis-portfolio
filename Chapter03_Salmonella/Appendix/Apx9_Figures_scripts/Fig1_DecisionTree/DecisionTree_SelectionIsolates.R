
setwd("~/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 03 Salmonella/WRITING_CHAPTER03/Figures_scripts/Fig_DecisionTree")


#install.packages("DiagrammeRsvg")

library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)

diagram <- grViz("
digraph nanopore_selection {

  graph [layout = dot, rankdir = LR, nodesep = 0.6, ranksep = 0.8]  # LR for side-by-side

  # Global node style
  node [shape=box, style=filled, fontname='Helvetica', fontsize=10]

  # Root
  start [label='Total sequencing capacity available: 451 isolates', fillcolor='#e0e0e0']

  # Two serovars
  bov [label='Bovismorbificans', fillcolor='#e3f2fd']
  give [label='Give', fillcolor='#d1c4e9']

  # Human isolates boxes
  bov_human [label='Selection rule: Human isolates\\n Select 1 per health district for each year', fillcolor='#e3f2fd']
  bov_human_out [label='Recovered: 121\\nNot-recovered: 10', shape=oval, fillcolor='#c8e6c9']

  give_human [label='Selection rule: Human isolates\\nSelect 1 per health district for each year', fillcolor='#d1c4e9']
  give_human_out [label='Recovered: 24\\nNot-recovered: 16', shape=oval, fillcolor='#c8e6c9']

  # Non-human isolates boxes (Bovismorbificans)
  bov_bovine [label='Selection rule: Bovine\\nIf >60/year reported → Select 25 per year\\nIf ≤60/year reported → Select 1 per year', fillcolor='#e3f2fd']
  bov_canine [label='Selection rule: Canine\\nIf ≥5/year reported → Select 2 per year\\nIf <5/year reported → Select 1 per year', fillcolor='#e3f2fd']
  bov_feline [label='Selection rule: Feline\\nIf ≥6/year reported → Select 2 per year\\nIf <6/year reported → Select 1 per year', fillcolor='#e3f2fd']
  bov_other [label='Selection rule: Other animals\\nSelect 1 per source per year', fillcolor='#e3f2fd']

  bov_nonhuman_out [label='Recovered: 234\\nNot-recovered: 16', shape=oval, fillcolor='#c8e6c9']

  # Non-human isolates boxes (Give)
  give_bovine [label='Selection rule: Bovine\\nIf <10/year reported → Select 1 per year\\nIf 10–50/year reported → Select 8 per year\\nIf ≥50/year reported → Select 10 per year', fillcolor='#d1c4e9']
  give_canine [label='Selection rule: Canine\\nIf 2–10/year reported → Select 3 per year\\nIf >10/year reported → Select 5 per year\\nIf <2/year reported → Select 1 per year', fillcolor='#d1c4e9']
  give_other [label='Selection rule: Other animals\\nSelect 1 per source per year', fillcolor='#d1c4e9']

  give_nonhuman_out [label='Recovered: 54\\nNot-recovered: 3', shape=oval, fillcolor='#c8e6c9']

  # Extras yellow box
  extras [label='433 recovered isolates: 18 remaining slots\\n Isolates added to complete flow cells\\n→ 6 extra Bovismorbificans (human)\\n→ 12 extra Bovismorbificans (bovine)', shape=box, fillcolor='#fff3e0']

  # Totals blue boxes
  bov_total [label='Sequenced: Bovismorbificans\\nHuman: 125\\nNon-human: 246\\nNot confirmed as Bovismorbificans: 4', fillcolor='#bbdefb']
  give_total [label='Sequenced: Give\\nHuman: 24\\nNon-human: 54\\nNot confirmed as Give: 1', fillcolor='#d1c4e9']
  
  # Final selection
  bov_analyses [label='Final selection for analyses: Bovismorbificans\\nTotal: 369', fillcolor='#bbdefb']
  give_analyses [label='Final selection for analyses: Give\\nTotal: 77', fillcolor='#d1c4e9']

  # Edges
  start -> bov
  start -> give

  # Human branches
  bov -> bov_human
  bov_human -> bov_human_out
  give -> give_human
  give_human -> give_human_out

  # Non-human branches (directly to green summary boxes)
  bov -> bov_bovine
  bov -> bov_canine
  bov -> bov_feline
  bov -> bov_other

  bov_bovine -> bov_nonhuman_out
  bov_canine -> bov_nonhuman_out
  bov_feline -> bov_nonhuman_out
  bov_other -> bov_nonhuman_out

  give -> give_bovine
  give -> give_canine
  give -> give_other

  give_bovine -> give_nonhuman_out
  give_canine -> give_nonhuman_out
  give_other -> give_nonhuman_out

  # Green boxes to yellow extras box
  bov_nonhuman_out -> extras
  bov_human_out -> extras
  give_nonhuman_out -> extras
  give_human_out -> extras

  # Yellow extras box to totals
  extras -> bov_total
  extras -> give_total
  
  #total analysed
  
  bov_total -> bov_analyses
  give_total -> give_analyses
  
  
}
")

diagram

# Convert to SVG
svg <- export_svg(diagram)

# Save as SVG
writeLines(svg, "decision_tree.svg")
