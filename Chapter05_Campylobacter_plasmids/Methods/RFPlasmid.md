# RFPlasmid Analysis

## Overview
This tool was used to predict plasmid and chromosome contigs from my two datasets:
D01: Assemblies generated with short-reads
D02: Polished assemblies generated with long-reads and short-reads
Only contigs predicted as plasmids with >60% probability (votes plasmid) were considered plasmid contigs.

GitHub page: <https://github.com/aldertzomer/RFPlasmid>


### Input Data
- Assembled genomes in FASTA format

#### RFPlasmid script

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=rfplasmid_2026 # job name
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=48:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=6 # number of cpus/threads
#SBATCH --mem=10G            # RAM

module purge
module load RFPlasmid/0.0.18-gimkl-2022a-Python-3.10.5

rfplasmid --species Campylobacter --input /nesi/project/massey03742/Campy_SACNZ/Illumina_assemblies/fasta --jelly --threads 6 --out /nesi/project/massey03742/Campy_SACNZ/RFPlasmid_Illumina

rfplasmid --species Campylobacter --input /nesi/project/massey03742/Campy_SACNZ/polished_assemblies/fasta --jelly --threads 6 --out /nesi/project/massey03742/Campy_SACNZ/RFPlasmid_Nanopore

```


---
