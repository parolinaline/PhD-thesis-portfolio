# Plasmid typing scheme

## Overview

Here I explain how I used the plasmid marker classification scheme developed by our collaborator to classify our plasmid contigs into different types based on the detection of these markers.
I used ABRICATE in order to detect the plasmid markers.

### Step-by-Step Execution
### Step 1: [Ruuning ABRICATE to detect markers]
#### Input Data
The input data will be the folder containing all the genome assemblies.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=abricate_plasmids_massey03742 # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=2G            # RAM
#SBATCH --cpus-per-task=4 # number of cpus/threads, stick to 2 for now
#SBATCH --time=10:00:00 # Walltime (HH:MM:SS)

cd /nesi/nobackup/massey03742/Campy_SACNZ/polished_assemblies/

module purge
module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1

#LOOP

for fasta_file in /nesi/nobackup/massey03742/Campy_SACNZ/polished_assemblies/*.fasta; do
    base_name=$(basename "$fasta_file" .fasta)
    abricate --datadir /nesi/nobackup/massey03742/plasmid_database/db/ --db CampyPlasmidtyping_abricate_20240710 "$fasta_file" > "${base_name}.tsv"
    echo "Finished processing $fasta_file"
done
```
#### Output

The output will be individual tsv files with the plasmid markers detected for each isolate. I used the following script to merge all tsv files together.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=tsv_together_massey03742 # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mem=1G            # RAM
#SBATCH --cpus-per-task=2 
#SBATCH --time=10:00:00 # Walltime (HH:MM:SS)

cd /nesi/nobackup/massey03742/Campy_SACNZ/polished_assemblies/

# Output file
output_file="combined_results_plasmid_campy_nanopore_all.txt"

# Initialize the output file
> "$output_file"

# Loop through all .tsv files
for tsv_file in *.tsv; do
    # Check if the file contains only a header or is empty
    if [[ $(wc -l < "$tsv_file") -le 1 ]]; then
        # If the file is empty or contains only the header, append a line to the output file
        echo -e "$tsv_file\tNo genes detected" >> "$output_file"
    else
        # If the file contains more than just the header, append its content to the output file
        cat "$tsv_file" >> "$output_file"
    fi
done
```

