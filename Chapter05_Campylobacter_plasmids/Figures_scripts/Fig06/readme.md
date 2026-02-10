## Description

Figure 6 – Horizontal visualisation and comparison of nucleotide sequences from conjugative plasmid types CP4 (A), CP5 (B), and CP10 (C). The cluster labels next to the isolate’s ID refer to plasmids with ≥80% of identity based on cd-hit analysis. The coloured arrows refer to the function of the genes according to Bakta annotation, as described in the legend. The grey links display the similarity of plasmid nucleotide sequences using BLASTn. GGGenomes   was used for visualisation.

### Input

1. All fasta and gff3 files from the plasmids that will be compared using GGGenomes.
2. Blastn results comparing the plasmid sequences in a tsv format
[Blastn_results](./plasmidM4_blast_filtered.tsv)  


## Blastn script
### Plasmid CP4

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_db_pM4
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/pM4

module purge
module load BLAST/2.16.0-GCC-12.3.0

#Make the database
makeblastdb -in plasmids_M4.fasta -dbtype nucl -out plasmid_M4_campy

# Run BLAST with multiple threads
blastn -query plasmids_M4.fasta \
       -db plasmid_M4_campy \
       -outfmt 6 \
       -out plasmidM4_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmidM4_blast.tsv > plasmidM4_blast_filtered.tsv

```
### Plasmid CP5

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_run_pM5
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/pM5

module purge
module load BLAST/2.16.0-GCC-12.3.0

#Make the database
makeblastdb -in plasmids_M5.fasta -dbtype nucl -out plasmid_M5_campy


# Run BLAST with multiple threads
blastn -query plasmids_M5.fasta \
       -db plasmid_M5_campy \
       -outfmt 6 \
       -out plasmidM5_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmidM5_blast.tsv > plasmidM5_blast_filtered.tsv
```
### Plasmid CP10

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_db_pM4
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/pM10

module purge
module load BLAST/2.16.0-GCC-12.3.0

#Make the database
makeblastdb -in plasmids_M10.fasta -dbtype nucl -out plasmid_M10_campy

# Run BLAST with multiple threads
blastn -query plasmids_M10.fasta \
       -db plasmid_M10_campy \
       -outfmt 6 \
       -out plasmidM10_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmidM10_blast.tsv > plasmidM10_blast_filtered.tsv
```
### BLASTn output
[Plasmids CP4 BLASTn](./plasmidM4_blast_filtered.tsv)  
[Plasmids CP5 BLASTn](./plasmidM5_blast_filtered.tsv)  
[Plasmids CP10 BLASTn](./plasmidM10_blast_filtered.tsv)

## GGGenomes
### Plasmids CP4 + CP5 + CP10

Now, we can use GGGenomes to visualise and annotate the plasmids comparison:

```r

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/pM4/")
getwd()

#load gggenomes and ggtree
library(gggenomes)
library(ggtree)
library(tidyverse)
library(randomcoloR)
library(aplot)

################################################################################
#==============================================================================#
#                         PLASMID GROUP CP4                                    #
#==============================================================================#
################################################################################

# Define paths to your files
plasmid_fasta_m4 <- c(
  "SC0014" = "./SC0014_39597_circlator.fasta",
  "SC0058" = "./SC0058_40676_circlator.fasta",
  "SC0084" = "./SC0084_38820_circlator.fasta",
  "SC0158" = "./SC0158_43614_circlator.fasta",
  "SC0212" = "./SC0212_42058_circlator.fasta",
  "SC0354" = "./SC0354_39851_circlator.fasta",
  "SC0391" = "./SC0391_36479_circlator.fasta",
  "SC0600" = "./SC0600_40695_circlator.fasta",
  "SC0947" = "./SC0947_46278_circlator.fasta",
  "SC0987" = "./SC0987_36277_circlator.fasta",
  "SC1210" = "./SC1210_42113_circlator.fasta",
  "SC1667" = "./SC1667_38016_circlator.fasta",
  "SC1704" = "./SC1704_40177_circlator.fasta",
  "SC1706" = "./SC1706_41836_circlator.fasta"
)

plasmid_gffs_m4 <- c(
  "SC0014" = "./SC0014B_39597_circlator.gff3",
  "SC0058" = "./SC0058_40676_circlator.gff3",
  "SC0084" = "./SC0084_38820_circlator.gff3",
  "SC0158" = "./SC0158_43614_circlator.gff3",
  "SC0212" = "./SC0212B_42058_circlator.gff3",
  "SC0354" = "./SC0354_39851_circlator.gff3",
  "SC0391" = "./SC0391_36479_circlator.gff3",
  "SC0600" = "./SC0600B_40695_circlator.gff3",
  "SC0947" = "./SC0947_46278_circlator.gff3",
  "SC0987" = "./SC0987C_36277_circlator.gff3",
  "SC1210" = "./SC1210B_42113_circlator.gff3",
  "SC1667" = "./SC1667C_38016_circlator.gff3",
  "SC1704" = "./SC1704_40177_circlator.gff3",
  "SC1706" = "./SC1706B_41836_circlator.gff3"
)

# Read sequences from FASTA files
seqs_m4 <- read_seqs(plasmid_fasta_m4)

# Check what was loaded
seqs_m4  

# From GFF3 files
genesM4 <- read_feats(plasmid_gffs_m4)

# Check the genes
genesM4

# Keep only protein-coding genes
genesM4 <- genesM4 %>%
  filter(type == "CDS")

# Check what columns you have
colnames(genesM4)

# Consolidate gene names - use 'gene' if available, otherwise use 'name'
genesM4 <- genesM4 %>%
  mutate(gene = coalesce(gene, name))

# Create the gggenomes object
pM4 <- gggenomes(
  genes = genesM4,
  seqs = seqs_m4
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqs_m4$seq_id)

# Check what seq_ids are in your genes
unique(genesM4$seq_id)


# Add layers
pM4 +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill = strand)) +
  geom_gene_tag(aes(label = name), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt6
linksM4 <- read_links("plasmidM4_blast_filtered.tsv", format = "blast")

head(linksM4)
colnames(linksM4)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_map_M4 <- c(
  "SC0014_39597" = "CP4_cluster0",
  "SC0058_40676" = "CP4_cluster1",
  "SC0084_38820" = "CP4_cluster2",
  "SC0158_43614" = "CP4_cluster3",
  "SC0212_42058" = "CP4_cluster4",
  "SC0354_39851" = "CP4_cluster5",
  "SC0391_36479" = "CP4_cluster6",
  "SC0600_40695" = "CP4_cluster7",
  "SC0947_46278" = "CP4_cluster8",
  "SC0987_36277" = "CP4_cluster9",
  "SC1210_42113" = "CP4_cluster10",
  "SC1667_38016" = "CP4_cluster11",
  "SC1704_40177" = "CP4_cluster12",
  "SC1706_41836" = "CP4_cluster13"
)

# Apply to seqs
seqs_m4 <- seqs_m4 %>%
  mutate(seq_id = label_map_M4[seq_id])

# Apply to genes
genesM4 <- genesM4 %>%
  mutate(seq_id = label_map_M4[seq_id])

# Apply to links
linksM4 <- linksM4 %>%
  mutate(
    seq_id = label_map_M4[seq_id],
    seq_id2 = label_map_M4[seq_id2]
  )


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genesM4 <- genesM4 %>%
  mutate(gene_category = case_when(
    gene %in% c("virB10","virB11","virD4","virB2","virB4","virB5","virB6","virB8","virB9") ~ "T4SS-like genes",
    gene == "mobC" ~ "Mobilisation",
    gene %in% c("ssb", "Replication gene") ~ "Replication",
    TRUE ~ "Others"
  ))

# Define colours
gene_coloursM4 <- c(
  "T4SS-like genes" = "steelblue",
  "Mobilisation" = "green",
  "Replication" = "orange",
  "Others" = "grey70"
)

# Recreate gggenomes object with updated genes
pM4 <- gggenomes(
  genes = genesM4,
  seqs = seqs_m4,
  links = linksM4
)


# Plot
final_plotM4 <- pM4 +
  geom_seq() +
  geom_bin_label(size = 4.5, hjust = 1) +
  geom_link(alpha = 0.3) +
  geom_gene(aes(fill = gene_category), size = 4) +
  geom_gene_tag(
    aes(label = ifelse(bin_id == first(bin_id), name, "")),
    nudge_y = 0.15,
    check_overlap = TRUE,
    size = 3.5,
    angle = 30
  ) +
  scale_fill_manual(
    values = gene_coloursM4,
    name = "Gene function",  # Legend title
    guide = "none"
  ) +
  ggtitle("A) CP4") +
  theme(plot.title = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

final_plotM4


################################################################################
#==============================================================================#
#                         PLASMID GROUP CP5                                    #
#==============================================================================#
################################################################################


getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/pM5/")
getwd()


# Define paths to your files
plasmid_fastasM5 <- c("SC0158"="./SC0158_30067_circlator.fasta",
                      "SC0592"="./SC0592_27689_circlator.fasta",
                      "SC0987"="./SC0987_28370_circlator.fasta",
                      "SC1586"="./SC1586_28955_circlator.fasta",
                      "SC1597"="./SC1597_27545_circlator.fasta",
                      "SC1667"="./SC1667_28032_circlator.fasta",
                      "SC0436"="./SC0436_26570_OG_circlator.fasta",
                      "SC1718"="./SC1718_25645_OG_circlator.fasta")

plasmid_gffsM5 <- c("SC0158"="./SC0158_30067_circlator.gff3",
                    "SC0592"="./SC0592A_27689_circlator.gff3",
                    "SC0987"="./SC0987A_28370_circlator.gff3",
                    "SC1586"="./SC1586_28955_circlator.gff3",
                    "SC1597"="./SC1597_27545_circlator.gff3",
                    "SC1667"="./SC1667A_28032_circlator.gff3",
                    "SC0436"="./SC0436_26570_OG_circlator.gff3",
                    "SC1718"="./SC1718_25645_OG_circlator.gff3")

# Read sequences from FASTA files
seqsM5 <- read_seqs(plasmid_fastasM5)

# Check what was loaded
seqsM5  

# From GFF3 files
genesM5 <- read_feats(plasmid_gffsM5)

# Check the genes
genesM5


# Create the gggenomes object
pM5 <- gggenomes(
  genes = genesM5,
  seqs = seqsM5
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqsM5$seq_id)

# Check what seq_ids are in your genes
unique(genesM5$seq_id)


# Add layers
pM5 +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill = strand)) +
  geom_gene_tag(aes(label = gene), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt6
linksM5 <- read_links("plasmidM5_blast_filtered.tsv", format = "blast")

head(linksM5)
colnames(linksM5)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_mapM5 <- c("SC0158_pM5"="CP5_cluster0",
                 "SC0436_pM5"="CP5_cluster6",
                 "SC0592_pM5"="CP5_cluster4",
                 "SC0987_pM5"="CP5_cluster2",
                 "SC1586_pM5"="CP5_cluster1",
                 "SC1597_pM5"="CP5_cluster5",
                 "SC1667_pM5"="CP5_cluster3",
                 "SC1718_pM5"="CP5_cluster7")

# Apply to seqs
seqsM5 <- seqsM5 %>%
  mutate(seq_id = label_mapM5[seq_id])

# Apply to genes
genesM5 <- genesM5 %>%
  mutate(seq_id = label_mapM5[seq_id])

# Apply to links
linksM5 <- linksM5 %>%
  mutate(
    seq_id = label_mapM5[seq_id],
    seq_id2 = label_mapM5[seq_id2]
  )


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genesM5 <- genesM5 %>%
  mutate(gene_category = case_when(
    gene %in% c("virB10","virB11","virD4","virB2","virB3","virB4","virB5","virB6","virB7","virB8","virB9") ~ "T4SS-like genes",
    gene == "mobA" ~ "Mobilisation",
    gene %in% c("ssb", "Rep", "Rep init", "repL") ~ "Replication",
    TRUE ~ "Others"
  ))

# Define colours
gene_coloursM5 <- c(
  "T4SS-like genes" = "steelblue",
  "Mobilisation" = "green",
  "Replication" = "orange",
  "Others" = "grey70"
)

# Recreate gggenomes object with updated genes
pM5 <- gggenomes(
  genes = genesM5,
  seqs = seqsM5,
  links = linksM5
)


# Plot

final_plot_M5 <- pM5 +
  geom_seq() +
  geom_bin_label(size = 4.5, hjust = 1) +
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
    values = gene_coloursM5,
    name = "Gene function",  # Legend title
    guide = "none"
  ) +
  ggtitle("B) CP5") +
  theme(plot.title = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

final_plot_M5

################################################################################
#==============================================================================#
#                         PLASMID GROUP CP4                                    #
#==============================================================================#
################################################################################

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/pM10/")
getwd()

# Define paths to your files
plasmid_fastasM10 <- c("SC0853_cluster0"="./SC0853_36918_circlator.fasta",
                       "SC0987_cluster3"="./SC0987_34827_circlator.fasta",
                       "SC1521_cluster2"="./SC1521_36038_circlator.fasta",
                       "SC1622_cluster1"="./SC1622_36560_circlator.fasta",
                       "SC1718_cluster4"="./SC1718_35240_circlator.fasta",
                       "SC1706_cluster5"="./SC1706_34734.fasta")

# Define paths to your files
plasmid_gffsM10 <- c("SC0853_cluster0"="./SC0853_36918_circlator.gff3",
                     "SC0987_cluster3"="./SC0987B_34827_circlator.gff3",
                     "SC1521_cluster2"="./SC1521_36038_circlator.gff3",
                     "SC1622_cluster1"="./SC1622B_36560_circlator.gff3",
                     "SC1718_cluster4"="./SC1718B_35240_circlator.gff3",
                     "SC1706_cluster5"="./SC1706_34734.gff3")

# Check unique strand values in your data
genes_rawM10 <- read_feats(plasmid_gffsM10)
unique(genes_rawM10$strand)

# Read sequences from FASTA files
seqsM10 <- read_seqs(plasmid_fastasM10)

# Check what was loaded
seqsM10  

# From GFF3 files
genesM10 <- read_feats(plasmid_gffsM10)

# Check the genes
genesM10


# Create the gggenomes object
pM10 <- gggenomes(
  genes = genesM10,
  seqs = seqsM10
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqsM10$seq_id)

# Check what seq_ids are in your genes
unique(genesM10$seq_id)


# Add layers
pM10 +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill = strand)) +
  geom_gene_tag(aes(label = gene), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt
linksM10 <- read_links("plasmidM10_blast_filtered.tsv", format = "blast")

head(linksM10)
colnames(linksM10)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_mapM10 <- c("SC0853_cluster0"="CP10_cluster0",
                  "SC0987_cluster3"="CP10_cluster3",
                  "SC1521_cluster2"="CP10_cluster2",
                  "SC1622_cluster1"="CP10_cluster1",
                  "SC1718_cluster4"="CP10_cluster4",
                  "SC1706_cluster5"="CP10_cluster5")

# Apply to seqs
seqsM10 <- seqsM10 %>%
  mutate(seq_id = label_mapM10[seq_id])

# Apply to genes
genesM10 <- genesM10 %>%
  mutate(seq_id = label_mapM10[seq_id])

# Apply to links
linksM10 <- linksM10 %>%
  mutate(
    seq_id = label_mapM10[seq_id],
    seq_id2 = label_mapM10[seq_id2]
  )


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genesM10 <- genesM10 %>%
  mutate(gene_category = case_when(
    gene %in% c("virB10","virB11","virD4","virB2","virB3","virB4","virB5","virB6","virB8","virB9","virB7") ~ "T4SS-like genes",
    gene %in% c("mobA", "mobC") ~ "Mobilisation",
    gene %in% c("ssb", "Replication gene", "repE") ~ "Replication",
    TRUE ~ "Others"
  ))

# Define colours
gene_coloursM10 <- c(
  "T4SS-like genes" = "steelblue",
  "Mobilisation" = "green",
  "Replication" = "orange",
  "Others" = "grey70"
)

# Recreate gggenomes object with updated genes
pM10 <- gggenomes(
  genes = genesM10,
  seqs = seqsM10,
  links = linksM10
)

# Plot
final_plot_M10 <- pM10 +
  geom_seq() +
  geom_bin_label(size = 4.5, hjust = 1) +
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
    values = gene_coloursM10,
    name = "Gene function"
  ) +
  ggtitle("C) CP10") +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "right")



final_plot_M10

################################################################################
###              PLOTING EVERYTHING TOGETHER                   #################
################################################################################

final_plot_ALL <- final_plotM4 %>%
  insert_bottom(final_plot_M5) %>%
  insert_bottom(final_plot_M10)

final_plot_ALL


# Save as PDF
ggsave("plasmid_ALL.pdf", plot = final_plot_ALL, width = 15, height = 20, units = "in")

```

### Final GGGenomes Figure

[GGGenomes_plot](./all_plasmids.jpg)
