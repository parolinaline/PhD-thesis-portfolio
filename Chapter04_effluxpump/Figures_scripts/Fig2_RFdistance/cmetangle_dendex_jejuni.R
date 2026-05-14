## ============================================================
## cme gene tree comparison vs cgMLST reference tree
## Campylobacter jejuni - efflux pump allelic diversity chapter
## ============================================================
## Workflow:
##   1. Load trees
##   2. Find shared tips, report excluded isolates
##   3. Prune all trees to shared tip set
##   4. Plot tanglegrams (cgMLST vs each cme gene tree)
##   5. Compute normalized Robinson-Foulds distances
## ============================================================

library(ape)
library(phangorn)
library(dendextend)
library(tidyverse)

setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 04 Campy cmeABC efflux pump/WRITING/Figures/tanglegram/")


## ── 1. Load trees ─────────────────────────────────────────────

snp_tree  <- read.tree("cgMLST_jejuni.newick")   # update filename if different
cmeA_tree <- read.tree("cmeA_jej.treefile")
cmeB_tree <- read.tree("cmeB_jej.treefile")
cmeC_tree <- read.tree("cmeC_jej.treefile")
cmeR_tree <- read.tree("cmeR_jej.treefile")

gene_trees <- list(
  cmeA = cmeA_tree,
  cmeB = cmeB_tree,
  cmeC = cmeC_tree,
  cmeR = cmeR_tree
)

## ── 2. Find shared tips & report excluded isolates ────────────

cat("\n=== TIP LABEL SUMMARY ===\n")
cat("cgMLST tree tips:", length(snp_tree$tip.label), "\n")
for (gene in names(gene_trees)) {
  cat(gene, "tree tips:", length(gene_trees[[gene]]$tip.label), "\n")
}

all_tip_sets <- c(
  list(SNP = snp_tree$tip.label),
  lapply(gene_trees, function(t) t$tip.label)
)

shared_tips <- Reduce(intersect, all_tip_sets)

cat("\n=== SHARED TIPS ACROSS ALL TREES ===\n")
cat("Isolates shared across all 5 trees:", length(shared_tips), "\n")

cat("\n=== EXCLUDED ISOLATES (present in cgMLST tree but missing from gene tree) ===\n")
for (gene in names(gene_trees)) {
  excluded <- setdiff(snp_tree$tip.label, gene_trees[[gene]]$tip.label)
  if (length(excluded) == 0) {
    cat(gene, ": no isolates excluded\n")
  } else {
    cat(gene, "(", length(excluded), "excluded ):\n")
    cat(paste(" -", excluded, collapse = "\n"), "\n")
  }
}

cat("\n=== ISOLATES IN GENE TREES BUT NOT IN cgMLST TREE (unexpected) ===\n")
for (gene in names(gene_trees)) {
  extra <- setdiff(gene_trees[[gene]]$tip.label, snp_tree$tip.label)
  if (length(extra) == 0) {
    cat(gene, ": none\n")
  } else {
    cat(gene, "(", length(extra), "unexpected ):\n")
    cat(paste(" -", extra, collapse = "\n"), "\n")
  }
}

## ── 3. Prune all trees to shared tip set ──────────────────────

prune_to_shared <- function(tree, keep_tips) {
  drop_tips <- setdiff(tree$tip.label, keep_tips)
  if (length(drop_tips) > 0) tree <- drop.tip(tree, drop_tips)
  return(tree)
}

snp_pruned  <- prune_to_shared(snp_tree, shared_tips)
gene_pruned <- lapply(gene_trees, prune_to_shared, keep_tips = shared_tips)

cat("\n=== AFTER PRUNING ===\n")
cat("All trees now have", length(shared_tips), "tips\n")

## ── 4. Convert phylo → dendrogram ─────────────────────────────
## dendextend works with dendrogram objects, not phylo directly.
## We convert via hclust using cophenetic distances, which preserves
## the tree topology for visualisation purposes.

phylo_to_dend <- function(tree) {
  tree %>%
    cophenetic() %>%          # pairwise distance matrix from tree
    as.dist() %>%
    hclust(method = "average") %>%
    as.dendrogram()
}

snp_dend  <- phylo_to_dend(snp_pruned)
gene_dend <- lapply(gene_pruned, phylo_to_dend)

## ── 5. Plot tanglegrams ────────────────────────────────────────
## One PDF per gene, cgMLST tree on the left, gene tree on the right.
## Connecting lines = identity (same isolate on both sides).
## Tip labels suppressed (too many isolates to read).

for (gene in names(gene_dend)) {

  outfile <- paste0("tanglegram_cgMLST_vs_", gene, "_cjejuni.pdf")
  pdf(outfile, width = 10, height = 14)

  # Only the RIGHT tree (gene tree) is rotated to minimise crossings.
  # After untangling, cgMLST is locked back to its original order
  # to guarantee identical layout across all four figures.
  dend_list <- dendlist(snp_dend, gene_dend[[gene]]) %>%
    untangle(method = "step1side")
  dend_list[[1]] <- snp_dend   # lock cgMLST back to original order

  tanglegram(
    dend_list,

    ## ── Labels ──────────────────────────────────────────────
    lab.cex        = 0.001,    # effectively invisible tip labels
    edge.lwd       = 1,        # branch line width

    ## ── Connecting lines ────────────────────────────────────
    lwd            = 0.6,      # connecting line width
    color_lines    = "grey50", # single colour for identity lines
    # tip labels are matched by name automatically

    ## ── Tip alignment lines ──────────────────────────────────
    dLeaf          = -0.01,    # negative value hides alignment dashes
    ## if dashes still visible, also set:
    # edge.root.lwd  = 0,

    ## ── Headers ─────────────────────────────────────────────
    main           = paste("cgMLST tree  vs ", gene, " tree — C. jejuni"),
    main_left      = "cgMLST",
    main_right     = gene,
    cex.main       = 1.1,

    ## ── General ─────────────────────────────────────────────
    margin_inner   = 5,        # space between the two trees for the lines
    margin_outer   = 1,
    axes           = FALSE     # suppress axes
  )

  dev.off()
  cat("Saved:", outfile, "\n")
}


cat("\nAll tanglegrams done.\n")

## ── 6. Robinson-Foulds distances (normalized) ─────────────────
## RF distance normalized by max possible RF = 2*(n_tips - 3) for unrooted
## Gives a value between 0 (identical topology) and 1 (maximally different)

cat("\n=== ROBINSON-FOULDS DISTANCES (normalized) ===\n")
cat(sprintf("%-10s %10s %10s %10s\n", "Gene", "RF_raw", "RF_max", "RF_norm"))

rf_results <- data.frame(
  gene    = character(),
  RF_raw  = integer(),
  RF_max  = integer(),
  RF_norm = numeric()
)

n_tips <- length(shared_tips)
RF_max <- 2 * (n_tips - 3)

for (gene in names(gene_pruned)) {

  rf_raw <- RF.dist(snp_pruned, gene_pruned[[gene]],
                    normalize = FALSE,
                    rooted    = FALSE)

  rf_norm <- rf_raw / RF_max

  cat(sprintf("%-10s %10d %10d %10.4f\n", gene, rf_raw, RF_max, rf_norm))

  rf_results <- rbind(rf_results, data.frame(
    gene    = gene,
    RF_raw  = rf_raw,
    RF_max  = RF_max,
    RF_norm = round(rf_norm, 4)
  ))
}

write.csv(rf_results, "cjejuni_RF_distances_cgMLST.csv", row.names = FALSE)
cat("\nRF distances saved to: cjejuni_RF_distances_cgMLST.csv\n")

## ── 7. Summary bar plot of RF distances ───────────────────────

rf_plot <- ggplot(rf_results, aes(x = gene, y = RF_norm, fill = gene)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = round(RF_norm, 3)), vjust = -0.4, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Normalized Robinson-Foulds distance vs cgMLST tree",
    subtitle = "C. jejuni — cme efflux pump genes",
    x        = "Gene tree",
    y        = "Normalized RF distance"
  ) +
  theme_classic(base_size = 13)

rf_plot
ggsave("cjejuni_RF_barplot_cgMLST.pdf", rf_plot, width = 7, height = 5)
cat("RF barplot saved to: cjejuni_RF_barplot_cgMLST.pdf\n")
