```bash 

#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=mob_suite_campy # job name
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=24:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=4 # number of cpus/threads
#SBATCH --mem=8G            # RAM

module purge
module load Apptainer/1.3.1

#This variable will assist with keeping the actual functional call short.

export CMD="apptainer exec /nesi/project/massey03742/software/mob_suite/mob_suite_3.1.9.aimg mob_typer"


#${CMD} --multi --infile #/nesi/nobackup/massey03742/Campy_Nanopore_SACNZ/best_assemblies/plasmid_contigs/plasmid_combined.fasta --out_file #Campy_mobsuite_typer_plasmids.txt 


${CMD} --multi --infile ./SC0134_2_np2.fasta --out_file Campy_mobtyper_SC0134.txt

```
