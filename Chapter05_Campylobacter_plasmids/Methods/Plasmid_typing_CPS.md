# Plasmid typing of circularisable plasmids using MashTree and Roary

## Overview

Here I explain the steps involved in classifying plasmid contigs generated with Illumina short reads dataset that were categorised as circularisable by SKESA.
I compare the whole-sequence of these circularisable plasmid sequences (CPS) using MashTree and Roary.
This refers to section 2.1.5 of my thesis.

## Workflow Summary

1. Get plasmid sequences that were predicted as plasmid contigs by RFPlasmid and circularisable by SKESA.
2. Extract these plasmid sequences and save in new individual fast files.
3. Use MashTree to compare these plasmid sequences and obtain the MashDistance.
4. Annotate these plasmids using Prokka.
5. Use gff3 files as input for Roary to obtain a distance matrix based on gene presence/absence data and calculate a distance matrix of these plasmid sequences using Jaccard.
6. Compare and stablish plasmid types based on both results.

### Step-by-Step Execution
### Step 1: [Selecting the plasmid contigs for the analysis]
#### Assemblies being used for this analysis
Only the assemblies from BioProject PRJNA675916 are being used for this analysis. These are assemblies generated with SKESA with Illumina short-reads.

#### Input Data
Use RFPlasmid output and SKESA output to select the contigs predicted as plasmids and categorised as circularisable by SKESA, which is an evidence of a complete sequence.
Generate a file that has the Contig ID and the name of the fasta file of the assembly where this contig is located. The whole-sequence file should be named after the isolate's ID.
Example:

INSERT EXAMPLE

### Step 2: [Extraction of selected contig sequences as individual fasta files]

```bash
