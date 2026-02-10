## Description

Figure 9 – Horizontal maps of representative pTet plasmids from each cd-hit cluster. The coloured arrows represent the gene's functional classification. The lines indicate the % identity of the amino acid sequences according to BLASTp. The plasmid comparison was generated and visualised using Clinker.

### Input

1. All fasta and gff3 files from the plasmids that will be compared using GGGenomes.
2. Blastn results comparing the plasmid sequences in a tsv format


## Blastn script
### Plasmid pTET SC1704 x NCBI nt database

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_db_pTET
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/Global_pTET

module load BLAST/2.16.0-GCC-12.3.0

makeblastdb -in plasmids_pTET.fasta -dbtype nucl -out plasmid_pTET_campy


# Run BLAST with multiple threads
blastn -query plasmids_pTET.fasta \
       -db plasmid_pTET_campy \
       -outfmt 6 \
       -out plasmid_pTET_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmid_pTET_blast.tsv > plasmid_pTET_blast_filt.tsv


```

### BLASTn output
[Plasmids pTET SC1704 BLASTn vs NCBI BLASTn database](./plasmid_pTET_blast_filt.tsv)  

## GGGenomes
### Plasmids pTET SC1704 BLASTn vs NCBI BLASTn database

Now, we can use GGGenomes to visualise and annotate the plasmids comparison:

```r

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/Global_pTET/")
getwd()

#load gggenomes and ggtree
library(gggenomes)
library(ggtree)
library(tidyverse)
library(randomcoloR)


# =============================================================================
# STEP 1: Define your file paths for pM5
# =============================================================================

plasmid_fastas <- c("SC1704_pTET"="./SC1704_PTET_circlator.fasta",
                  "AP025988_pTET"="./AP025988.1_circlator.fasta",
                  "CP030889_pTET"="./CP030889.1_circlator.fasta",
                  "AP028406_pTET"="./AP028406.1_circlator.fasta",
                  "AY394560_pTET"="./AY394560.1_circlator.fasta",
                  "CP011017_pTET"="./CP011017.1_circlator.fasta",
                  "CP017861_pTET"="./CP017861.1_circlator.fasta",
                  "CP048762_pTET"="./CP048762.1_circlator.fasta",
                  "CP172365_pTET"="./CP172365.1_circlator.fasta",
                  "ON101335_pTET"="./ON101335.1_circlator.fasta",
                  "OQ553941_pTET"="./OQ553941.1_circlator.fasta"
                  )

plasmid_gffs <- c("SC1704_pTET"="./SC1704_tetO.gff3",
                "AP025988_pTET"="./AP025988.1_circlator.gff3",
                "CP030889_pTET"="./CP030889.1_circlator.gff3",
                "AP028406_pTET"="./AP028406.1_circlator.gff3",
                "AY394560_pTET"="./AY394560.1_circlator.gff3",
                "CP011017_pTET"="./CP011017.1_circlator.gff3",
                "CP017861_pTET"="./CP017861.1_circlator.gff3",
                "CP048762_pTET"="./CP048762.1_circlator.gff3",
                "CP172365_pTET"="./CP172365.1_circlator.gff3",
                "ON101335_pTET"="./ON101335.1_circlator.gff3",
                "OQ553941_pTET"="./OQ553941.1_circlator.gff3"
                )

# Check unique strand values in your data
genes_raw <- read_feats(plasmid_gffs)
unique(genes_raw$strand)

# Read sequences from FASTA files
seqs <- read_seqs(plasmid_fastas)

# Check what was loaded
seqs  

# From GFF3 files
genes <- read_feats(plasmid_gffs)

# Check the genes
genes


# Create the gggenomes object
p <- gggenomes(
genes = genes,
seqs = seqs
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqs$seq_id)

# Check what seq_ids are in your genes
unique(genes$seq_id)


# Add layers
p +
geom_seq() +
geom_bin_label() +
geom_gene(aes(fill = strand)) +
geom_gene_tag(aes(label = gene), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt
links <- read_links("plasmid_pTET_blast_filt.tsv", format = "blast")

head(links)
colnames(links)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_map <- c("AP025988_pTET"="AP025988_pTET",
             "CP030889_pTET"="CP030889_pTET",
             "AP028406_pTET"="AP028406_pTET",
             "AY394560_pTET"="AY394560_pTET",
             "CP011017_pTET"="CP011017_pTET",
             "CP017861_pTET"="CP017861_pTET",
             "CP048762_pTET"="CP048762_pTET",
             "CP172365_pTET"="CP172365_pTET",
             "ON101335_pTET"="ON101335_pTET",
             "OQ553941_pTET"="OQ553941_pTET",
             "SC1704_pTET"="SC1704_pTET")

# Apply to seqs
seqs <- seqs %>%
mutate(seq_id = label_map[seq_id])

# Apply to genes
genes <- genes %>%
mutate(seq_id = label_map[seq_id])

# Apply to links
links <- links %>%
mutate(
  seq_id = label_map[seq_id],
  seq_id2 = label_map[seq_id2]
)


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genes <- genes %>%
mutate(gene_category = case_when(
  gene %in% c("virB10","virB11","virD4","virB2","virB3","virB4","virB5","virB6","virB7","virB8","virB9") ~ "T4SS-like genes",
  gene  %in% c("mobC","mobA", "mob") ~ "Mobilisation",
  gene %in% c("ssb", "Replication gene", "rep") ~ "Replication",
  gene %in% c("tetO","PTet37","PTet38","tetL", "tet(O)") ~ "Tetracycline resistance",
  TRUE ~ "Others"
))

# Define colours
gene_colours <- c(
"T4SS-like genes" = "steelblue",
"Mobilisation" = "green",
"Replication" = "orange",
"Tetracycline resistance" = "darkred",
"Others" = "grey70"
)

# Recreate gggenomes object with updated genes
p <- gggenomes(
genes = genes,
seqs = seqs,
links = links
)

# Plot
final_plot_pTET <- p +
geom_seq() +
geom_bin_label(size = 5, hjust = 1) +
geom_link(alpha = 0.3) +
geom_gene(aes(fill = gene_category), size = 4) +
geom_gene_tag(
  aes(label = ifelse(bin_id == first(bin_id), gene, "")),
  nudge_y = 0.15,
  check_overlap = TRUE,
  size = 3.5,
  angle = 30
) +
scale_fill_manual(
  values = gene_colours,
  name = "Gene function"  # Legend title
) +
theme(legend.position = "bottom")

final_plot_pTET


# Save as PDF
ggsave("plasmid_pTET_plot.pdf", plot = final_plot_pTET, width = 13, height = 10, units = "in")


```

### Final GGGenomes Figure

[GGGenomes_plot](./pTET_GGGenomes.jpg)
