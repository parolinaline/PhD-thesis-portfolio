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


${CMD} search --pointfinder-organism campylobacter -o starAMR_campy /nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/polished_assemblies/*.fasta

```

