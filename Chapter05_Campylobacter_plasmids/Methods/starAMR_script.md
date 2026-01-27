# starAMR pipeline

## Overview
This tool was used to detect point mutations associated with antimicrobial resistance phenotype (pointfinder), resistance genes using ResFinder, and confirm the ST using the MLST database.

GitHub page: <(https://github.com/phac-nml/staramr)>


### Input Data
- Assembled genomes in FASTA format

#### starAMR script

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=StarAMR_massey03742
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=120:00:00

module purge
module load Apptainer/1.3.1

#This variable will assist with keeping the actual functional call short. Change it from "bakta" to other commands
#such as bakta_db as needed
export CMD="apptainer exec /nesi/project/massey03742/software/containers/staramr-0.10.0.aimg staramr"


${CMD} search --pointfinder-organism campylobacter -o starAMR_campy /nesi/nobackup/massey03742/Campy_SACNZ/polished_assemblies/*.fasta


```

---
