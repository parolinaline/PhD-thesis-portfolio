## Description

Figure 8 - Horizontal visualisation and comparison of nucleotide sequences from conjugative plasmid types CP4 (A), CP5 (B), and CP10 (C) retrieved from the NCBI database nt. and compared to representative isolates from this study. The coloured arrows refer to the function of the genes according to Bakta annotation, as described in the legend. The grey links display the similarity of plasmid nucleotide sequences using BLASTn. GGGenomes was used for visualisation.

### Input

1. All fasta and gff3 files from the plasmids that will be compared using GGGenomes.
2. Blastn results comparing the plasmid sequences in a tsv format


## Blastn script
### Plasmid CP4 x NCBI nt database

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_db_GpM4
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/Global_pM4

module load BLAST/2.16.0-GCC-12.3.0

makeblastdb -in plasmids_GM4.fasta -dbtype nucl -out plasmid_GM4_campy



# Run BLAST with multiple threads
blastn -query plasmids_GM4.fasta \
       -db plasmid_GM4_campy \
       -outfmt 6 \
       -out plasmid_GM4_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmid_GM4_blast.tsv > plasmid_GM4_blast_filt.tsv


```
### Plasmid CP5 x NCBI nt database

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_db_GpM5
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/Global_pM5

module purge
module load BLAST/2.16.0-GCC-12.3.0


makeblastdb -in plasmids_GpM5.fasta -dbtype nucl -out plasmid_GpM5_campy


# Run BLAST with multiple threads
blastn -query plasmids_GpM5.fasta \
       -db plasmid_GpM5_campy \
       -outfmt 6 \
       -out plasmid_GpM5_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmid_GpM5_blast.tsv > plasmid_GpM5_blast_filt.tsv

```
### Plasmid CP10 x NCBI nt database

```bash
#!/bin/bash -e
#SBATCH --job-name=blast_db_GpM5
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

cd /nesi/project/massey03742/Campy_Nanopore_SACNZ/BLASTn_gggenomes/Global_pM10

module purge
module load BLAST/2.16.0-GCC-12.3.0


makeblastdb -in plasmids_GpM10.fasta -dbtype nucl -out plasmid_GpM10_campy



# Run BLAST with multiple threads
blastn -query plasmids_GpM10.fasta \
       -db plasmid_GpM10_campy \
       -outfmt 6 \
       -out plasmid_GpM10_blast.tsv \
       -evalue 1e-10 \
       -num_threads 4

# Filter self-hits
awk '$1 != $2' plasmid_GpM10_blast.tsv > plasmid_GpM10_blast_filt.tsv

```
### BLASTn output
[Plasmids CP4 BLASTn vs NCBI BLASTn database](./plasmid_GM4_blast_filt.tsv)  
[Plasmids CP5 BLASTn vs NCBI BLASTn database](./plasmid_GpM5_blast_filt.tsv)  
[Plasmids CP10 BLASTn vs NCBI BLASTn database](./plasmid_GpM10_blast_filt.tsv)

## GGGenomes
### Plasmids CP4 + CP5 + CP10 vs NCBI BLASTn database

Now, we can use GGGenomes to visualise and annotate the plasmids comparison:

```r


#load gggenomes and ggtree
library(gggenomes)
library(ggtree)
library(tidyverse)
library(randomcoloR)
library(aplot)

################################################################################
#==============================================================================#
#                         PLASMID GROUP CP4 GLOBAL                             #
#==============================================================================#
################################################################################


getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/Global_pM4/")
getwd()


# =============================================================================
# STEP 1: Define your file paths for pM5
# =============================================================================

plasmid_fastasG4 <- c("SC0947_pM04"="./SC0947_46278_circlator.fasta",
                      "AP028347_pM04"="./AP028347.1_circlator.fasta",
                      "AP028368_pM04"="./AP028368.1_circlator.fasta",
                      "AY394560_pM04"="./AY394560.1_circlator.fasta",
                      "CP001961_pM04"="./CP001961.1_circlator.fasta",
                      "CP010073_pM04"="./CP010073.1_circlator.fasta",
                      "CP017419_pM04"="./CP017419.1_circlator.fasta",
                      "CP043764_pM04"="./CP043764.1_circlator.fasta",
                      "CP048762_pM04"="./CP048762.1_circlator.fasta",
                      "CP131443_pM04"="./CP131443.1_circlator.fasta",
                      "CP172365_pM04"="./CP172365.1_circlator.fasta")

plasmid_gffsG4 <- c("SC0947_pM04"="./SC0947_46278_circlator.gff3",
                    "AP028347_pM04"="./AP028347.1_circlator.gff3",
                    "AP028368_pM04"="./AP028368.1_circlator.gff3",
                    "AY394560_pM04"="./AY394560.1_circlator.gff3",
                    "CP001961_pM04"="./CP001961.1_circlator.gff3",
                    "CP010073_pM04"="./CP010073.1_circlator.gff3",
                    "CP017419_pM04"="./CP017419.1_circlator.gff3",
                    "CP043764_pM04"="./CP043764.1_circlator.gff3",
                    "CP048762_pM04"="./CP048762.1_circlator.gff3",
                    "CP131443_pM04"="./CP131443.1_circlator.gff3",
                    "CP172365_pM04"="./CP172365.1_circlator.gff3")

# Check unique strand values in your data
genes_rawG4 <- read_feats(plasmid_gffsG4)
unique(genes_rawG4$strand)

# Read sequences from FASTA files
seqsG4 <- read_seqs(plasmid_fastasG4)

# Check what was loaded
seqsG4  

# From GFF3 files
genesG4 <- read_feats(plasmid_gffsG4)

# Check the genes
genesG4


# Create the gggenomes object
pG4 <- gggenomes(
  genes = genesG4,
  seqs = seqsG4
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqsG4$seq_id)

# Check what seq_ids are in your genes
unique(genesG4$seq_id)


# Add layers
pG4 +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill = strand)) +
  geom_gene_tag(aes(label = gene), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt
linksG4 <- read_links("plasmid_GM4_blast_filt.tsv", format = "blast")

head(linksG4)
colnames(linksG4)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_mapG4 <- c("SC0947_pM04"="SC0947_pM04",
                 "AP028347_pM04"="AP028347_pM04",
                 "AP028368_pM04"="AP028368_pM04",
                 "AY394560_pM04"="AY394560_pM04",
                 "CP001961_pM04"="CP001961_pM04",
                 "CP010073_pM04"="CP010073_pM04",
                 "CP017419_pM04"="CP017419_pM04",
                 "CP043764_pM04"="CP043764_pM04",
                 "CP048762_pM04"="CP048762_pM04",
                 "CP131443_pM04"="CP131443_pM04",
                 "CP172365_pM04"="CP172365_pM04")

# Apply to seqs
seqsG4 <- seqsG4 %>%
  mutate(seq_id = label_mapG4[seq_id])

# Apply to genes
genesG4 <- genesG4 %>%
  mutate(seq_id = label_mapG4[seq_id])

# Apply to links
linksG4 <- linksG4 %>%
  mutate(
    seq_id = label_mapG4[seq_id],
    seq_id2 = label_mapG4[seq_id2]
  )


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genesG4 <- genesG4 %>%
  mutate(gene_category = case_when(
    gene %in% c("virB10","virB11","virD4","virB2","virB3","virB4","virB5","virB6","virB7","virB8","virB9") ~ "T4SS-like genes",
    gene  %in% c("mobC","mobA", "mob") ~ "Mobilisation",
    gene %in% c("ssb", "Replication gene", "rep") ~ "Replication",
    gene %in% c("tetO","PTet37","PTet38","tetL") ~ "Tetracycline resistance",
    TRUE ~ "Others"
  ))

# Define colours
gene_coloursG4 <- c(
  "T4SS-like genes" = "steelblue",
  "Mobilisation" = "green",
  "Replication" = "orange",
  "Tetracycline resistance" = "darkred",
  "Others" = "grey70"
)

# Recreate gggenomes object with updated genes
pG4 <- gggenomes(
  genes = genesG4,
  seqs = seqsG4,
  links = linksG4
)


# Plot
final_plot_GPM4 <- pG4 +
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
    values = gene_coloursG4,
    name = "Gene function",  # Legend title
    guide = "none"
  ) +
  ggtitle("A) Global comparison of CP4 plasmids") +
  theme(plot.title = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

final_plot_GPM4

################################################################################
#==============================================================================#
#                         PLASMID GROUP CP5 GLOBAL                             #
#==============================================================================#
################################################################################


getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/Global_pM5/")
getwd()



# =============================================================================
# STEP 1: Define your file paths for pM5
# =============================================================================

plasmid_fastasG5 <- c("SC0600_pM10"="./SC0600_29581_M5_circlator.fasta",
                      "CP044172_pM10"="./CP044172.1_circlator.fasta",
                      "CP066747_pM10"="./CP066747.1_circlator.fasta",
                      "CP071588_pM10"="./CP071588.1_circlator.fasta",
                      "CP092035_pM10"="./CP092035.1_circlator.fasta",
                      "CP125385_pM10"="./CP125385.1_circlator.fasta",
                      "CP132117_pM10"="./CP132117.1_circlator.fasta",
                      "CP132125_pM10"="./CP132125.1_circlator.fasta",
                      "CP132127_pM10"="./CP132127.1_circlator.fasta",
                      "CP132129_pM10"="./CP132129.1_circlator.fasta",
                      "CP132131_pM10"="./CP132131.1_circlator.fasta")


plasmid_gffsG5 <- c("SC0600_pM10"="./SC0600_29581_M5_circlator.gff3",
                    "CP044172_pM10"="./CP044172.1_circlator.gff3",
                    "CP066747_pM10"="./CP066747.1_circlator.gff3",
                    "CP071588_pM10"="./CP071588.1_circlator.gff3",
                    "CP092035_pM10"="./CP092035.1_circlator.gff3",
                    "CP125385_pM10"="./CP125385.1_circlator.gff3",
                    "CP132117_pM10"="./CP132117.1_circlator.gff3",
                    "CP132125_pM10"="./CP132125.1_circlator.gff3",
                    "CP132127_pM10"="./CP132127.1_circlator.gff3",
                    "CP132129_pM10"="./CP132129.1_circlator.gff3",
                    "CP132131_pM10"="./CP132131.1_circlator.gff3")

# Check unique strand values in your data
genes_rawG5 <- read_feats(plasmid_gffsG5)
unique(genes_rawG5$strand)

# Read sequences from FASTA files
seqsG5 <- read_seqs(plasmid_fastasG5)

# Check what was loaded
seqsG5  

# From GFF3 files
genesG5 <- read_feats(plasmid_gffsG5)

# Check the genes
genesG5


# Create the gggenomes object
pG5 <- gggenomes(
  genes = genesG5,
  seqs = seqsG5
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqsG5$seq_id)

# Check what seq_ids are in your genes
unique(genesG5$seq_id)


# Add layers
pG5 +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill = strand)) +
  geom_gene_tag(aes(label = gene), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt
linksG5 <- read_links("plasmid_GpM5_blast_filt.tsv", format = "blast")

head(linksG5)
colnames(linksG5)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_mapG5 <- c("SC0600_pM10"="SC0600_pM05",
                 "CP044172_pM10"="CP044172_pM05",
                 "CP066747_pM10"="CP066747_pM05",
                 "CP071588_pM10"="CP071588_pM05",
                 "CP092035_pM10"="CP092035_pM05",
                 "CP125385_pM10"="CP125385_pM05",
                 "CP132117_pM10"="CP132117_pM05",
                 "CP132125_pM10"="CP132125_pM05",
                 "CP132127_pM10"="CP132127_pM05",
                 "CP132129_pM10"="CP132129_pM05",
                 "CP132131_pM10"="CP132131_pM05")

# Apply to seqs
seqsG5 <- seqsG5 %>%
  mutate(seq_id = label_mapG5[seq_id])

# Apply to genes
genesG5 <- genesG5 %>%
  mutate(seq_id = label_mapG5[seq_id])

# Apply to links
linksG5 <- linksG5 %>%
  mutate(
    seq_id = label_mapG5[seq_id],
    seq_id2 = label_mapG5[seq_id2]
  )


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genesG5 <- genesG5 %>%
  mutate(gene_category = case_when(
    gene %in% c("virB10","virB11","virD4","virB2","virB3","virB4","virB5","virB6","virB7","virB8","virB9") ~ "T4SS-like genes",
    gene == "mobA" ~ "Mobilisation",
    gene %in% c("ssb", "Replication gene", "repL") ~ "Replication",
    TRUE ~ "Others"
  ))

# Define colours
gene_coloursG5 <- c(
  "T4SS-like genes" = "steelblue",
  "Mobilisation" = "green",
  "Replication" = "orange",
  "Others" = "grey70"
)

# Recreate gggenomes object with updated genes
pG5 <- gggenomes(
  genes = genesG5,
  seqs = seqsG5,
  links = linksG5
)

# Plot
final_plot_GPM5 <- pG5 +
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
    values = gene_coloursG5,
    name = "Gene function",  # Legend title
    guide = "none"
  ) +
  ggtitle("B) Global comparison of CP5 plasmids") +
  theme(plot.title = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

final_plot_GPM5

################################################################################
#==============================================================================#
#                         PLASMID GROUP CP10 GLOBAL                            #
#==============================================================================#
################################################################################

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 05 Campy plasmids/gggenomes/Global_M10/")
getwd()

# =============================================================================
# STEP 1: Define your file paths for pM5
# =============================================================================

plasmid_fastasG10 <- c("SC0853_pM10"="./SC0853_pM10_circlator.fasta",
                       "CP010074_pM10"="./CP010074.1_circlator.fasta",
                       "CP038864_pM10"="./CP038864.1_circlator.fasta",
                       "CP116956_pM10"="./CP116956.1_circlator.fasta",
                       "CP140347_pM10"="./CP140347.1_circlator.fasta",
                       "CP140351_pM10"="./CP140351.1_circlator.fasta",
                       "CP140363_pM10"="./CP140363.1_circlator.fasta",
                       "CP140387_pM10"="./CP140387.1_circlator.fasta",
                       "CP140399_pM10"="./CP140399.1_circlator.fasta",
                       "CP181391_pM10"="./CP181391.1_circlator.fasta",
                       "OQ553945_pM10"="./OQ553945.1_circlator.fasta")

plasmid_gffsG10 <- c("SC0853_pM10"="./SC0853_pM10_circlator.gff3",
                     "CP010074_pM10"="./CP010074.1_circlator.gff3",
                     "CP038864_pM10"="./CP038864.1_circlator.gff3",
                     "CP116956_pM10"="./CP116956.1_circlator.gff3",
                     "CP140347_pM10"="./CP140347.1_circlator.gff3",
                     "CP140351_pM10"="./CP140351.1_circlator.gff3",
                     "CP140363_pM10"="./CP140363.1_circlator.gff3",
                     "CP140387_pM10"="./CP140387.1_circlator.gff3",
                     "CP140399_pM10"="./CP140399.1_circlator.gff3",
                     "CP181391_pM10"="./CP181391.1_circlator.gff3",
                     "OQ553945_pM10"="./OQ553945.1_circlator.gff3")

# Check unique strand values in your data
genes_rawG10 <- read_feats(plasmid_gffsG10)
unique(genes_rawG10$strand)

# Read sequences from FASTA files
seqsG10 <- read_seqs(plasmid_fastasG10)

# Check what was loaded
seqsG10  

# From GFF3 files
genesG10 <- read_feats(plasmid_gffsG10)

# Check the genes
genesG10


# Create the gggenomes object
pG10 <- gggenomes(
  genes = genesG10,
  seqs = seqsG10
)

#let's check if the headers match, they have to!
# Check what seq_ids are in your sequences
unique(seqsG10$seq_id)

# Check what seq_ids are in your genes
unique(genesG10$seq_id)


# Add layers
pG10 +
  geom_seq() +
  geom_bin_label() +
  geom_gene(aes(fill = strand)) +
  geom_gene_tag(aes(label = gene), nudge_y = 0.1, check_overlap = TRUE)


# Read the BLAST links
# Specify the format as BLAST outfmt
linksG10 <- read_links("plasmid_GpM10_blast_filt.tsv", format = "blast")

head(linksG10)
colnames(linksG10)

################################################################################
###                 NEW LABELS FOR EACH SEQUENCE                       #########
################################################################################

# Create the mapping
label_mapG10 <- c("SC0853_pM10"="SC0853_pM10",
                  "CP010074_pM10"="CP010074_pM10",
                  "CP038864_pM10"="CP038864_pM10",
                  "CP116956_pM10"="CP116956_pM10",
                  "CP140347_pM10"="CP140347_pM10",
                  "CP140351_pM10"="CP140351_pM10",
                  "CP140363_pM10"="CP140363_pM10",
                  "CP140387_pM10"="CP140387_pM10",
                  "CP140399_pM10"="CP140399_pM10",
                  "CP181391_pM10"="CP181391_pM10",
                  "OQ553945_pM10"="OQ553945_pM10")

# Apply to seqs
seqsG10 <- seqsG10 %>%
  mutate(seq_id = label_mapG10[seq_id])

# Apply to genes
genesG10 <- genesG10 %>%
  mutate(seq_id = label_mapG10[seq_id])

# Apply to links
linksG10 <- linksG10 %>%
  mutate(
    seq_id = label_mapG10[seq_id],
    seq_id2 = label_mapG10[seq_id2]
  )


################################################################################
###               COLOURING GENES                              #################
################################################################################


# Create gene categories
genesG10 <- genesG10 %>%
  mutate(gene_category = case_when(
    gene %in% c("virB10","virB11","virD4","virB2","virB3","virB4","virB5","virB6","virB7","virB8","virB9", "T4SS", "trbC", "conjugative transfer protein") ~ "T4SS-like genes",
    gene == "mobA" ~ "Mobilisation",
    gene %in% c("ssb", "Replication gene", "repE") ~ "Replication",
    TRUE ~ "Others"
  ))

# Define colours
gene_coloursG10 <- c(
  "T4SS-like genes" = "steelblue",
  "Mobilisation" = "green",
  "Replication" = "orange",
  "Others" = "grey70"
)

# Recreate gggenomes object with updated genes
pG10 <- gggenomes(
  genes = genesG10,
  seqs = seqsG10,
  links = linksG10
)


# Plot
final_plot_GPM10 <- pG10 +
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
    values = gene_coloursG10,
    name = "Gene function",  # Legend title
  ) +
  ggtitle("C) Global comparison of CP10 plasmids") +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "right")


final_plot_GPM10


################################################################################
###              PLOTING EVERYTHING TOGETHER                   #################
################################################################################

final_plot_ALL_GLOBAL <- final_plot_GPM4 %>%
  insert_bottom(final_plot_GPM5) %>%
  insert_bottom(final_plot_GPM10)

final_plot_ALL_GLOBAL


# Save as PDF
ggsave("plasmid_GLOBAL_ALL.pdf", plot = final_plot_ALL_GLOBAL, width = 15, height = 25, units = "in")


```

### Final GGGenomes Figure

[GGGenomes_plot](./GGGenomes_Global.jpg)
