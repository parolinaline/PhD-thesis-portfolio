# cgMLST

## Overview
Here I show the step-by-step proccess to run the cgMLST and obtain a Mininum spanning tree based on the cgMLST output.
I am using the cgMLST schema for Campylobacter available on pubMLST.
I am using chewBACCA to run this analysis.
The fasta sequences used as input are the assemblies generated with Nanopore long reads and Illumina short reads. 

## Workflow Summary
1. Download the Campylobacter cgMLST schema from PubMLST (species ID 6, schema 1).
2. Prepares/validates an external schema (the one from pubMLST) for use with chewBBACA (converts FASTA files to proper format).
3. Run allele calling on your assemblies against the prepared schema.
4. Extracts the core genes (present in ≥95% of isolates) from the allele call results.


### Step-by-Step Execution
### Step 1: [Download the pubMLST schema for Campylobacter]
#### Script

```bash
#!/bin/bash -e
#SBATCH --job-name=chewbbaca_campy
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

module purge
module load chewBBACA/3.1.2-gimkl-2022a-Python-3.11.3

chewBBACA.py DownloadSchema -sp 6 -sc 1 -o /nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/best_assemblies/cgmlst_campy --latest
```
### Step 2: [Prepare the pubMLST schema to run with chewBBACA]
#### Script

```bash

#!/bin/bash -e
#SBATCH --job-name=chewbbaca_prep_schema
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

module purge
module load chewBBACA/3.1.2-gimkl-2022a-Python-3.11.3


chewBBACA.py PrepExternalSchema -i campylobacter_cgmlst_fastas/ -o campylobacter_pubMLST_prepared

```

### Step 3: [Run allele calling on assemblies]
#### Script

```bash
#!/bin/bash -e
#SBATCH --job-name=chewbbaca_run_campy
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

module purge
module load chewBBACA/3.1.2-gimkl-2022a-Python-3.11.3


chewBBACA.py AlleleCall \
  -i /nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/best_assemblies/fasta \
  -g /nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/best_assemblies/campylobacter_pubMLST_prepared/ \
  -o /nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/best_assemblies/pubMLST_results --ptf Campylobacter_jejuni.trn --cpu 8
```

### Step 4: [Extracts the core genes]
#### Script

```bash

#!/bin/bash -e
#SBATCH --job-name=chewbbaca_extract_cgMLST
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      

module purge
module load chewBBACA/3.1.2-gimkl-2022a-Python-3.11.3

chewBBACA.py ExtractCgMLST -i /nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/best_assemblies/cgMLST_pubMLST_results/results_alleles_NoParalogs.tsv \
  -o cgMLST_profiles_pubMLST \
  --t 0.95 
```

#### Output

Use cgMLST95.tsv as an input to generate a minimum spanning tree with GrapeTree
---
