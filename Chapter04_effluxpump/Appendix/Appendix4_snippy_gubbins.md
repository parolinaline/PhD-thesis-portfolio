# SNP Calling, Recombination Detection, and Visualisation of the *cmeABCR* Locus

## Overview

This page documents the bioinformatics workflow used to identify single nucleotide polymorphisms (SNPs), detect recombination events, and visualise recombination within the *cmeABC* efflux pump operon and its regulator *cmeR* in *Campylobacter jejuni* and *Campylobacter coli*. Short-read Illumina sequences were obtained from BioProject [PRJNA675916](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA675916). Analyses were run separately for *C. jejuni* and *C. coli* because species-specific reference genomes were used — complete hybrid assemblies (Nanopore + Illumina) generated from isolates within this project.

All computationally intensive steps were run on the NeSI high-performance computing platform using SLURM job scheduling. Only representative scripts are shown; equivalent scripts were run for both species.

---

## Workflow Summary

| Step | Tool | Purpose |
|------|------|---------|
| 1 | FastQC + MultiQC | Raw read quality assessment |
| 2 | Trimmomatic | Adapter trimming and quality filtering |
| 3 | FastQC + MultiQC | Post-trimming quality check |
| 4 | Snippy (via Nullarbor) | Per-isolate variant calling against reference |
| 5 | snippy-core | Generate whole-genome core SNP alignment |
| 6 | Gubbins | Recombination detection and masking |
| 7 | IQ-TREE | Maximum-likelihood phylogenetic inference |
| 8 | R | Extract recombination events overlapping *cmeABCR* |
| 9 | R | Visualise recombination aligned to phylogenetic tree |

---

## Step 1 — Raw Read Quality Assessment

Raw reads were assessed with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and summarised with [MultiQC](https://multiqc.info/). The script below loops over all paired-end FASTQ files, runs FastQC, collects basic per-sample statistics (read counts, average length, proportion of reads <50 bp, and file integrity), and writes a summary CSV.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=qc_all_reads
#SBATCH --output=slurmlog/qc_all_reads.%j.out
#SBATCH --error=slurmlog/qc_all_reads.%j.err
#SBATCH --time=96:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16

module purge
module load FastQC/0.12.1
module load MultiQC/1.24.1-foss-2023a-Python-3.11.6

FASTQ_DIR="/nesi/nobackup/massey03742/effluxpump_snippy"
FASTQC_OUT="${FASTQ_DIR}/fastqc_results"
STATS_OUT="${FASTQ_DIR}/read_stats"
MULTIQC_OUT="${FASTQ_DIR}/multiqc_report"

mkdir -p ${FASTQC_OUT} ${STATS_OUT} ${MULTIQC_OUT}

SUMMARY_CSV="${MULTIQC_OUT}/basic_stats_summary.csv"
echo "Sample,R1_size_MB,R2_size_MB,R1_reads,R2_reads,R1_avg_len,R2_avg_len,R1_short%,R2_short%,Status" > ${SUMMARY_CSV}

TOTAL=$(ls ${FASTQ_DIR}/SRR*_1.fastq.gz | wc -l)
COUNTER=0

echo "Starting QC analysis for ${TOTAL} samples"

for R1 in ${FASTQ_DIR}/SRR*_1.fastq.gz; do
    COUNTER=$((COUNTER + 1))
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE=$(basename ${R1} _1.fastq.gz)

    echo "Processing sample ${COUNTER}/${TOTAL}: ${SAMPLE}"

    fastqc -t ${SLURM_CPUS_PER_TASK} -o ${FASTQC_OUT} ${R1} ${R2}

    STATS_FILE="${STATS_OUT}/${SAMPLE}_stats.txt"
    R1_SIZE=$(du -m ${R1} | cut -f1)
    R2_SIZE=$(du -m ${R2} | cut -f1)
    R1_READS=$(zcat ${R1} | echo $(($(wc -l)/4)))
    R2_READS=$(zcat ${R2} | echo $(($(wc -l)/4)))

    if [ ${R1_READS} -eq ${R2_READS} ]; then STATUS="MATCHED"; else STATUS="MISMATCHED"; fi

    R1_AVG=$(zcat ${R1} | awk 'NR%4==2 {sum+=length($0); count++} END {printf "%.1f", sum/count}')
    R2_AVG=$(zcat ${R2} | awk 'NR%4==2 {sum+=length($0); count++} END {printf "%.1f", sum/count}')
    SHORT_R1=$(zcat ${R1} | awk 'NR%4==2 {if(length($0)<50) short++; total++} END {printf "%.2f", (short/total)*100}')
    SHORT_R2=$(zcat ${R2} | awk 'NR%4==2 {if(length($0)<50) short++; total++} END {printf "%.2f", (short/total)*100}')

    echo "${SAMPLE},${R1_SIZE},${R2_SIZE},${R1_READS},${R2_READS},${R1_AVG},${R2_AVG},${SHORT_R1},${SHORT_R2},${STATUS}" >> ${SUMMARY_CSV}
done

multiqc ${FASTQC_OUT} -o ${MULTIQC_OUT} -n all_samples_qc_report --force

echo "QC complete. MultiQC report: ${MULTIQC_OUT}/all_samples_qc_report.html"
```

---

## Step 2 — Read Trimming with Trimmomatic

Reads were trimmed using [Trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic) in paired-end mode. Nextera PE adapters were clipped, leading and trailing low-quality bases (Q < 3) were removed, a sliding window filter (4 bp window, minimum mean quality Q20) was applied, and reads shorter than 50 bp after trimming were discarded.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=trimmomatic_jejuni
#SBATCH --output=slurmlog/trimmomatic.%j.out
#SBATCH --error=slurmlog/trimmomatic.%j.err
#SBATCH --time=96:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

module purge
module load Trimmomatic/0.39-Java-1.8.0_144

RAW_DIR="/nesi/nobackup/massey03742/effluxpump_snippy"
TRIMMED_DIR="/nesi/nobackup/massey03742/effluxpump_snippy/jejuni_trimmed_reads"
ADAPTER_FILE="/nesi/nobackup/massey03742/effluxpump_snippy/adapters/NexteraPE-PE.fa"
THREADS=4

mkdir -p ${TRIMMED_DIR}/unpaired

TOTAL=$(ls ${RAW_DIR}/*_1.fastq.gz | wc -l)
COUNT=0

for R1 in ${RAW_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename ${R1} _1.fastq.gz)
    R2=${RAW_DIR}/${SAMPLE}_2.fastq.gz
    COUNT=$((COUNT + 1))

    # Skip already-processed samples
    if [ -f "${TRIMMED_DIR}/${SAMPLE}_1_paired.fastq.gz" ] && \
       [ -f "${TRIMMED_DIR}/${SAMPLE}_2_paired.fastq.gz" ]; then
        echo "[${COUNT}/${TOTAL}] ${SAMPLE} already processed, skipping."
        continue
    fi

    echo "[${COUNT}/${TOTAL}] Trimming ${SAMPLE}..."

    trimmomatic PE -threads ${THREADS} \
        ${R1} ${R2} \
        ${TRIMMED_DIR}/${SAMPLE}_1_paired.fastq.gz \
        ${TRIMMED_DIR}/unpaired/${SAMPLE}_1_unpaired.fastq.gz \
        ${TRIMMED_DIR}/${SAMPLE}_2_paired.fastq.gz \
        ${TRIMMED_DIR}/unpaired/${SAMPLE}_2_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:50

    echo "[${COUNT}/${TOTAL}] ${SAMPLE} done."
done

echo "All ${TOTAL} samples trimmed."
```

> **Post-trimming QC:** Read quality was re-assessed using the same FastQC/MultiQC script from Step 1 to confirm that trimming achieved the desired quality thresholds before proceeding.

---

## Step 3 — SNP Calling with Snippy

SNPs were called per isolate using [Snippy v4.6](https://github.com/tseemann/snippy), accessed via the Nullarbor module on NeSI. Separate reference genomes were used for *C. jejuni* (`SC0017_jejuni_ref.gbff`) and *C. coli* (`SC0001_coli_ref.gbff`) — both complete hybrid assemblies from this project.

### Step 3a — Generate the Snippy run script

`snippy-multi` takes a tab-separated input file (`ID`, `R1`, `R2`) and generates the per-isolate Snippy commands. The script below produces `runme_snippy_jej01.sh`:

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=01snippy_jejuni
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=5G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

module purge
module load nullarbor/2.0.20191013

snippy-multi input_fixed01.tab --ref SC0017_jejuni_ref.gbff --cpus 16 > runme_snippy_jej01.sh
```

The input tab file format is:

```
# ID    R1                        R2
SC0324  /path/to/SRR_1.fastq.gz   /path/to/SRR_2.fastq.gz
SC0326  /path/to/SRR_1.fastq.gz   /path/to/SRR_2.fastq.gz
```

### Step 3b — Run Snippy

The generated script (`runme_snippy_jej01.sh`) was opened with `nano` and SLURM headers and the Nullarbor module load were prepended before submission. The edited script runs one `snippy` call per isolate:

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=02jej_snippy
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00

module purge
module load nullarbor/2.0.20191013

snippy --outdir 'SC0324' \
    --R1 '/nesi/nobackup/massey03742/effluxpump_snippy/SRR17900495_1.fastq.gz' \
    --R2 '/nesi/nobackup/massey03742/effluxpump_snippy/SRR17900495_2.fastq.gz' \
    --ref jejuni_reference.gbff --cpus 16

snippy --outdir 'SC0326' \
    --R1 '/nesi/nobackup/massey03742/effluxpump_snippy/SRR17900493_1.fastq.gz' \
    --R2 '/nesi/nobackup/massey03742/effluxpump_snippy/SRR17900493_2.fastq.gz' \
    --ref jejuni_reference.gbff --cpus 16

# ... (continues for all isolates)

# Generate core SNP alignment from all per-isolate Snippy output directories
snippy-core --ref jejuni_reference.gbff SC0324 SC0326 # ... all isolate dirs
```

> `snippy-core` produces `clean.full.aln` (whole-genome alignment with invariant sites, used as Gubbins input) and `core.aln` (core SNP sites only). Both *C. jejuni* and *C. coli* datasets were processed identically.

---

## Step 4 — Recombination Detection with Gubbins

[Gubbins v3.2.2](https://github.com/nickjcroucher/gubbins) was used to identify and mask putative recombination regions in the core genome alignment. The `--remove-identical-sequences` flag was included to reduce computational overhead. After Gubbins finished, recombinant positions in the alignment were masked (converted to `N`) using `mask_gubbins_aln.py`, producing a recombination-free alignment for phylogenetic inference.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=jejuni_gubbins
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00

cd /nesi/nobackup/massey03742/effluxpump_snippy/snippy_jejuni

module purge
module load Gubbins/3.2.2-gimkl-2022a-Python-3.10.5
module load IQ-TREE/2.2.2.2-gimpi-2022a

INPUT_ALN=clean.full.aln
PREFIX=gubbins_out

# Run Gubbins — detects recombination using RAxML internally
run_gubbins.py -p ${PREFIX} ${INPUT_ALN} \
    --threads 8 \
    --remove-identical-sequences

echo "=== Gubbins complete ==="

# Mask recombinant positions in the alignment
mask_gubbins_aln.py \
    --aln ${INPUT_ALN} \
    --gff ${PREFIX}.recombination_predictions.gff \
    --out ${PREFIX}.masked.fasta

echo "=== Recombination masking complete ==="
```

**Key Gubbins outputs used downstream:**

| File | Description |
|------|-------------|
| `gubbins_out.recombination_predictions.gff` | Coordinates of predicted recombination blocks per taxon |
| `gubbins_out.final_tree.tre` | Phylogenetic tree (used for visualisation) |
| `gubbins_out.masked.fasta` | Recombination-masked alignment (used for IQ-TREE) |

---

## Step 5 — Phylogenetic Inference with IQ-TREE

A maximum-likelihood phylogenetic tree was inferred from the recombination-masked alignment using [IQ-TREE v2.2.2](http://www.iqtree.org/). IQ-TREE is bundled in the same module loaded in Step 4 and was run as part of the same job, or as a separate submission after Gubbins completed.

```bash
# Run after Gubbins masking — executed in the same environment
iqtree2 -s gubbins_out.masked.fasta \
    -m GTR+G \
    -B 1000 \
    -T 8 \
    --prefix iqtree_masked
```

> The GTR+G substitution model was used, with 1,000 ultrafast bootstrap replicates. The resulting tree (`iqtree_masked.treefile`) was used for all downstream visualisations.

---

## Step 6 — Extracting Recombination Events Overlapping *cmeABCR*

An R script was used to parse the Gubbins GFF output and identify recombination events that overlap the coordinates of the *cmeABC* operon and its regulator *cmeR*. Coordinates were taken from the species-specific reference genomes used in Steps 3–4 (verified with `grep -P "cme[ABCR]" ref.gff | awk '{print $4,$5,$7}'`). The script produces a per-isolate CSV table and a summary CSV indicating which genes were affected in each isolate.

This was run separately for *C. coli* and *C. jejuni* using their respective GFF files and gene coordinates.

### *C. coli*

Gene coordinates from reference SC0001:

| Gene | Start | End |
|------|-------|-----|
| *cmeC* | 373,168 | 374,646 |
| *cmeB* | 374,639 | 377,761 |
| *cmeA* | 377,761 | 378,864 |
| *cmeR* | 378,972 | 379,607 |

```r
library(tidyverse)

rec_gff_file <- "coli_gubbins_recombination.gff"
output_csv   <- "cmeABCRcoli_recombination_table.csv"

# Gene coordinates (C. coli reference SC0001)
gene_df <- tibble::tribble(
  ~gene,   ~gene_start, ~gene_end,
  "cmeC",  373168,      374646,
  "cmeB",  374639,      377761,
  "cmeA",  377761,      378864,
  "cmeR",  378972,      379607
)

# Parse Gubbins GFF — strip header lines and extract taxa per event
raw <- readLines(rec_gff_file)
raw <- raw[!grepl("^##", raw)]

rec_df <- read.table(text = raw, sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE, fill = TRUE, quote = "") %>%
  setNames(c("seq", "prog", "type", "start", "end",
             "score", "strand", "frame", "attributes")) %>%
  rowwise() %>%
  mutate(
    taxa_str = str_extract(attributes, "taxa=\\s*[^;]+"),
    taxa_str = str_remove(taxa_str, "taxa=\\s*"),
    taxa_str = str_trim(taxa_str),
    taxa     = list(
      str_split(taxa_str, "\\s+")[[1]] %>%
        str_remove_all('"') %>%
        .[nchar(.) > 0]
    ),
    n_taxa   = length(taxa),
    rec_type = if_else(n_taxa == 1, "unique", "shared")
  ) %>%
  ungroup() %>%
  select(start, end, taxa, n_taxa, rec_type)

# Filter recombinations overlapping the cmeABCR region
region_start <- min(gene_df$gene_start)
region_end   <- max(gene_df$gene_end)

rec_in_region <- rec_df %>%
  filter(end >= region_start & start <= region_end)

# Determine which genes each event overlaps
rec_with_genes <- rec_in_region %>%
  rowwise() %>%
  mutate(
    affected_genes = list(
      gene_df %>%
        filter(gene_end >= start & gene_start <= end) %>%
        pull(gene)
    ),
    genes_hit = paste(affected_genes, collapse = ", "),
    n_genes   = length(affected_genes)
  ) %>%
  ungroup()

# Unnest to one row per taxon per recombination event
result_table <- rec_with_genes %>%
  unnest(taxa) %>%
  select(
    isolate       = taxa,
    rec_start     = start,
    rec_end       = end,
    rec_type,
    n_taxa_in_rec = n_taxa,
    genes_hit,
    n_genes
  ) %>%
  mutate(rec_length = rec_end - rec_start + 1) %>%
  arrange(isolate, rec_start)

# Summary: one row per isolate
summary_table <- result_table %>%
  group_by(isolate) %>%
  summarise(
    n_recombinations = n(),
    total_rec_bp     = sum(rec_length),
    has_cmeC         = any(str_detect(genes_hit, "cmeC")),
    has_cmeB         = any(str_detect(genes_hit, "cmeB")),
    has_cmeA         = any(str_detect(genes_hit, "cmeA")),
    has_cmeR         = any(str_detect(genes_hit, "cmeR")),
    all_genes_hit    = paste(unique(unlist(str_split(genes_hit, ", "))),
                             collapse = ", "),
    rec_types        = paste(unique(rec_type), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_recombinations))

write.csv(result_table,  file = output_csv, row.names = FALSE)
write.csv(summary_table, file = gsub(".csv", "_summary.csv", output_csv),
          row.names = FALSE)

cat("\n=== cmeABCR Recombination Summary (C. coli) ===\n")
cat("Total recombination events in region:", nrow(rec_in_region), "\n")
cat("Isolates affected:", n_distinct(result_table$isolate), "\n")
cat("  cmeC:", sum(summary_table$has_cmeC), "isolates\n")
cat("  cmeB:", sum(summary_table$has_cmeB), "isolates\n")
cat("  cmeA:", sum(summary_table$has_cmeA), "isolates\n")
cat("  cmeR:", sum(summary_table$has_cmeR), "isolates\n")
```

### *C. jejuni*

Gene coordinates from reference SC0017:

| Gene | Start | End |
|------|-------|-----|
| *cmeC* | 342,811 | 344,289 |
| *cmeB* | 344,282 | 347,404 |
| *cmeA* | 347,404 | 348,507 |
| *cmeR* | 348,601 | 349,233 |

The same script was run with the *C. jejuni* GFF file (`gubbins_jejuni_predictions.gff`) and the coordinates above.

---

## Step 7 — Visualising Recombination in the *cmeABCR* Region

A multi-panel figure was produced in R combining: (i) a maximum-likelihood phylogenetic tree, (ii) a recombination segment plot aligned to the *cmeABCR* genomic region, (iii) gene annotation arrows, (iv) a per-position recombination density heatmap, and (v) isolate metadata tiles (ST, source). The figure was assembled using `aplot` and `cowplot`.

The full script is shown here for *C. coli*; an equivalent script was run for *C. jejuni* with species-specific input files and gene coordinates.

```r
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(gggenes)
library(aplot)
library(patchwork)
library(cowplot)
library(RColorBrewer)

# ── INPUT FILES ───────────────────────────────────────────────────────────────
tree_file    <- "./coli_gubbins_tree.tre"
rec_gff_file <- "coli_gubbins_recombination.gff"
meta_file    <- "Metadata_coli.csv"
output_file  <- "c_coli_cmeABCR_recombination.pdf"

# ── GENE COORDINATES (C. coli reference SC0001) ───────────────────────────────
gene_df <- tibble::tribble(
  ~gene,   ~start,  ~end,    ~strand, ~direction,
  "cmeC",  373168,  374646,  -1,      -1,
  "cmeB",  374639,  377761,  -1,      -1,
  "cmeA",  377761,  378864,  -1,      -1,
  "cmeR",  378972,  379607,  -1,      -1
)

padding      <- 500
region_start <- min(gene_df$start) - padding
region_end   <- max(gene_df$end)   + padding

output_width  <- 16
output_height <- 10

# ── COLOUR PALETTES ───────────────────────────────────────────────────────────
rec_colours <- c("shared" = "#D73027", "unique" = "#4575B4")
gene_fills  <- c("cmeR" = "#FDB462", "cmeA" = "#80B1D3",
                 "cmeB" = "#FB8072", "cmeC" = "#B3DE69")

# For large numbers of metadata categories (>8), interpolate across hues
many_colour_palette <- colorRampPalette(
  c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#999999",
    "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
    "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
)

###############################################################################
## 1. READ AND PROCESS TREE
###############################################################################

tree   <- ape::read.tree(tree_file)
n_taxa <- length(tree$tip.label)

# Auto-scale line and label sizes based on dataset size
tree_linewidth <- case_when(
  n_taxa < 20  ~ 0.8,
  n_taxa < 50  ~ 0.5,
  n_taxa < 100 ~ 0.35,
  TRUE         ~ 0.25
)

tip_label_size <- case_when(
  n_taxa < 20  ~ 3.5,
  n_taxa < 50  ~ 2.5,
  n_taxa < 100 ~ 2.0,
  TRUE         ~ 1.5
)

p_tree <- ggtree(tree, size = tree_linewidth) +
  geom_tiplab(size = tip_label_size, align = FALSE, offset = 0.0001,
              fontface = "plain") +
  theme_tree2() +
  theme(
    legend.position = "none",
    axis.line.x     = element_line(linewidth = 0.25),
    plot.margin     = margin(5, 0, 5, 5)
  )

taxa_order <- rev(get_taxa_name(p_tree))

###############################################################################
## 2. PARSE GUBBINS GFF
###############################################################################

parse_gubbins_gff <- function(gff_file) {
  raw <- readLines(gff_file)
  raw <- raw[!grepl("^##", raw)]

  df <- read.table(text = raw, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE, fill = TRUE, quote = "")
  colnames(df) <- c("seq", "prog", "type", "start", "end",
                    "score", "strand", "frame", "attributes")

  df %>%
    rowwise() %>%
    mutate(
      taxa_str = str_extract(attributes, "taxa=\\s*[^;]+"),
      taxa_str = str_remove(taxa_str, "taxa=\\s*"),
      taxa_str = str_trim(taxa_str),
      taxa     = list(str_split(taxa_str, "\\s+")[[1]])
    ) %>%
    ungroup() %>%
    mutate(
      n_taxa   = map_int(taxa, length),
      rec_type = if_else(n_taxa == 1, "unique", "shared")
    ) %>%
    select(start, end, taxa, n_taxa, rec_type)
}

rec_df <- parse_gubbins_gff(rec_gff_file)

# Clip recombination events to the plotting region
rec_region <- rec_df %>%
  filter(end >= region_start & start <= region_end) %>%
  mutate(
    start = pmax(start, region_start),
    end   = pmin(end, region_end)
  )

rec_long <- rec_region %>%
  unnest(taxa) %>%
  rename(label = taxa) %>%
  mutate(label = factor(label, levels = taxa_order))

rec_linewidth <- case_when(
  n_taxa < 20  ~ 5,
  n_taxa < 50  ~ 2.5,
  n_taxa < 100 ~ 1.5,
  TRUE         ~ 0.8
)

###############################################################################
## 3. RECOMBINATION PANEL
###############################################################################

p_rec <- ggplot(rec_long,
                aes(x = start, xend = end,
                    y = label, yend = label,
                    colour = rec_type)) +
  geom_segment(alpha = 0.6, linewidth = rec_linewidth) +
  scale_colour_manual(
    values = rec_colours,
    name   = "Recombination",
    labels = c("shared" = "Shared (≥2 taxa)", "unique" = "Unique (1 taxon)")
  ) +
  scale_x_continuous(
    limits = c(region_start, region_end),
    expand = c(0, 0),
    labels = scales::comma
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Genomic position (bp)") +
  theme_bw() +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    plot.margin      = margin(5, 10, 5, 0)
  )

###############################################################################
## 4. GENE ANNOTATION ARROWS (cmeABCR)
###############################################################################

gene_plot_df <- gene_df %>%
  mutate(molecule = "cmeABCR_operon", midpoint = (start + end) / 2)

p_genes <- ggplot(gene_plot_df,
                  aes(xmin = start, xmax = end, y = molecule,
                      fill = gene, forward = direction == 1)) +
  geom_gene_arrow(
    arrowhead_height  = unit(4, "mm"),
    arrowhead_width   = unit(2, "mm"),
    arrow_body_height = unit(3, "mm")
  ) +
  geom_gene_label(aes(label = gene), align = "centre", grow = FALSE,
                  height = grid::unit(4, "mm"), fontface = "italic") +
  scale_fill_manual(values = gene_fills, guide = "none") +
  scale_x_continuous(limits = c(region_start, region_end), expand = c(0, 0)) +
  theme_genes() +
  theme(
    axis.text   = element_blank(),
    axis.title  = element_blank(),
    axis.ticks  = element_blank(),
    axis.line   = element_blank(),
    panel.grid  = element_blank(),
    plot.margin = margin(2, 10, 0, 0)
  )

###############################################################################
## 5. RECOMBINATION DENSITY HEATMAP
###############################################################################

density_df <- rec_region %>%
  mutate(bases = map2(start, end, seq)) %>%
  select(bases) %>%
  unnest(bases) %>%
  count(bases, name = "count") %>%
  mutate(
    count = pmin(count, quantile(count, 0.95)),  # cap at 95th percentile
    y     = 1
  )

p_heatmap <- ggplot(density_df, aes(x = bases, y = y,
                                    fill = count, colour = count)) +
  geom_tile() +
  scale_x_continuous(limits = c(region_start, region_end), expand = c(0, 0)) +
  scale_fill_gradient2(
    low      = "navy", mid = "orange", high = "red",
    midpoint = max(density_df$count) / 2,
    name     = "Recombination\nevents",
    aesthetics = c("fill", "colour")
  ) +
  guides(
    fill   = guide_colorbar(title.position = "left", title.hjust = 0.5,
                            direction = "horizontal", barwidth = 8),
    colour = "none"
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "navy"),
    legend.position  = "bottom",
    plot.margin      = margin(0, 10, 0, 0)
  )

###############################################################################
## 6. METADATA PANEL (ST, Source)
###############################################################################

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
if ("Id" %in% names(meta)) meta <- rename(meta, id = Id)
if ("ID" %in% names(meta)) meta <- rename(meta, id = ID)

meta_cols <- c("ST", "Source")
meta_cols <- meta_cols[meta_cols %in% names(meta)]

if (length(meta_cols) == 0) {
  warning("No matching metadata columns found. Skipping metadata panel.")
  p_meta <- NULL

} else {

  meta_long <- meta %>%
    select(id, all_of(meta_cols)) %>%
    filter(id %in% taxa_order) %>%
    mutate(across(all_of(meta_cols), as.character)) %>%
    pivot_longer(-id, names_to = "variable", values_to = "value") %>%
    mutate(
      id       = factor(id, levels = taxa_order),
      variable = factor(variable, levels = meta_cols),
      value    = as.character(value)
    )

  palette_names <- c("Set2", "Set3", "Paired", "Dark2", "Accent")
  meta_palettes <- list()

  for (i in seq_along(meta_cols)) {
    col    <- meta_cols[i]
    vals   <- unique(meta_long$value[meta_long$variable == col])
    n_vals <- length(vals)

    if (n_vals > 8) {
      colours <- many_colour_palette(n_vals)
    } else {
      pal_name <- palette_names[((i - 1) %% length(palette_names)) + 1]
      colours  <- brewer.pal(max(3, n_vals), pal_name)[seq_len(n_vals)]
    }

    names(colours)       <- vals
    meta_palettes[[col]] <- colours
  }

  make_meta_tile <- function(col_name, pal) {
    df_sub <- meta_long %>% filter(variable == col_name)
    ggplot(df_sub, aes(x = variable, y = id, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = pal, name = col_name, na.value = "grey80") +
      scale_y_discrete(drop = FALSE) +
      theme_void() +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 8,
                                       face = "bold"),
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 8),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        plot.margin     = margin(5, 1, 5, 1)
      ) +
      guides(fill = guide_legend(
        title.position = "top",
        ncol           = min(4, length(pal)),
        title.hjust    = 0.5
      ))
  }

  meta_plots           <- map2(meta_cols, meta_palettes[meta_cols], make_meta_tile)
  meta_plots_no_legend <- map(meta_plots, ~ .x + theme(legend.position = "none"))
  meta_legends         <- map(meta_plots, ~ cowplot::get_legend(.x))
  p_meta               <- wrap_plots(meta_plots_no_legend, nrow = 1)
}

###############################################################################
## 7. ASSEMBLE FINAL FIGURE
###############################################################################

options("aplot_guides" = "keep")

if (!is.null(p_meta)) {
  combined <- p_rec %>%
    aplot::insert_left(p_meta, width = 0.15 * length(meta_cols)) %>%
    aplot::insert_left(p_tree, width = 0.6)
} else {
  combined <- p_rec %>%
    aplot::insert_left(p_tree, width = 0.6)
}

combined <- combined %>%
  aplot::insert_top(p_genes,  height = 0.06) %>%
  aplot::insert_top(p_heatmap + theme(legend.position = "none"), height = 0.03)

# Collect and arrange legends
rec_legend     <- cowplot::get_legend(p_rec     + theme(legend.position = "bottom"))
heatmap_legend <- cowplot::get_legend(p_heatmap + theme(legend.position = "bottom"))

if (!is.null(p_meta)) {
  meta_legends_explicit <- map2(meta_cols, meta_palettes[meta_cols], function(col_name, pal) {
    df_sub <- meta_long %>% filter(variable == col_name)
    p_leg <- ggplot(df_sub, aes(x = variable, y = id, fill = value)) +
      geom_tile() +
      scale_fill_manual(
        values = pal, name = col_name, na.value = "grey80",
        guide  = guide_legend(title.position = "top", title.hjust = 0.5,
                              ncol = min(5, length(pal)),
                              override.aes = list(size = 3))
      ) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 9),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.key      = element_rect(colour = "white", linewidth = 0.3)
      )
    cowplot::get_legend(p_leg)
  })
  all_legends <- c(list(rec_legend, heatmap_legend), meta_legends_explicit)
} else {
  all_legends <- list(rec_legend, heatmap_legend)
}

all_legends <- Filter(Negate(is.null), all_legends)

legend_row <- cowplot::plot_grid(
  plotlist   = all_legends,
  nrow       = 1,
  rel_widths = c(0.15, 0.2,
                 rep(0.65 / max(length(all_legends) - 2, 1),
                     max(length(all_legends) - 2, 0)))
)

final_plot <- cowplot::plot_grid(
  aplot::as.patchwork(combined),
  legend_row,
  nrow        = 2,
  rel_heights = c(1, 0.18)
)

###############################################################################
## 8. SAVE OUTPUT
###############################################################################

ggsave(
  filename = output_file,
  plot     = final_plot,
  width    = output_width,
  height   = output_height,
  dpi      = 300
)

cat("\n✓ Plot saved to:", output_file, "\n")
cat("  Region:", scales::comma(region_start), "-",
    scales::comma(region_end), "bp\n")
cat("  Taxa:", n_taxa, "\n")
cat("  Recombination events in region:", nrow(rec_region), "\n")
```

---

## Notes

1. **Gene coordinates** must be verified against your reference GFF before running:
   ```bash
   grep -P "cme[ABCR]" your_reference.gff | awk '{print $4,$5,$7}'
   ```

2. **Species-specific references:** *C. jejuni* analyses used reference SC0017 (NCTC 11168 locus tags: *cmeR* = Cj0368c, *cmeA* = Cj0367c, *cmeB* = Cj0366c, *cmeC* = Cj0365c). *C. coli* analyses used reference SC0001.

3. **Tip label scaling:** The R visualisation script auto-scales tip label and line sizes based on dataset size. Adjust `tip_label_size` or the tree panel `width` parameter in Section 7 if labels are clipped or overlapping.

4. **Colour palettes:** Metadata panels auto-scale: datasets with ≤8 categories use RColorBrewer palettes; datasets with >8 categories use a 16-hue interpolated palette (`many_colour_palette`). Missing or blank values are shown in grey.

5. **Metadata CSV format expected:**
   ```
   id,ST,Source
   SC0001,ST-827,Poultry
   SC0002,ST-828,Human
   ```
