# Assemblying genomes

## Overview

Here I explain how the assemblies were obtained for Chapter 05, starting from the filtering process of the reads until the polishing steps. I used two different approaches, and the final assemblies were compared and chosen as described on the thesis. The tools and versions used are indicated in the scripts. 

## Workflow Summary

1. Filtering Nanopore long-reads using FiltLong
2. Assemblying genomes:
    1. Unicycler hybrid assembly mode using Nanopore long-reads and Illumina short-reads
    2. Flye assembly using only Nanopore long-reads
3. Polish assemblies using Nanopore long-reads and Illumina short-reads 


### Step-by-Step Execution
#### Step 1: [Filtering Nanopore long-reads after basecalling]
#### Input Data

| File | Format | Description |
|------|--------|-------------|
| input_file.fastq | fastq | [Nanopore long-reads] |

I filtered the Nanopore long-reads using Filtlong.
<https://github.com/rrwick/Filtlong>
I used the following parameters:
```bash
--min_length,    I used 1000: Discard any read which is shorter than 1kbp
--keep_percent,    I used 95: It removes 5% of the worst reads. This is measured by bp, not by read count.
```

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=filtlong # job name (shows up in the queue)
#SBATCH --output=slurmlog/%A.%a.%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=10:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=6 # number of cpus/threads, stick to 2 for now
#SBATCH --mem=8G            # RAM

#access the working directory
cd /nesi/nobackup/massey03742/Campy_SACNZ/reads

module purge
module load Filtlong/0.2.0

# Loop through each porechop fastq file in the current directory
for input_file in *.fastq; do
    # Get the base name of the file (without the .fastq extension)
    base_name=$(basename "$input_file".fastq)
    # Define the output file name
    output_file="${base_name}_FL.fastq.gz"
    # Run Filtlong on the input file and save the output
    filtlong --min_length 1000 --keep_percent 95 "$input_file" | gzip > "$output_file"
done

```

#### Step 2: [Assemblying genomes]

I had the Illumina short-reads available for these isolates: BioProject PRJNA675916
I tested two assembly tools and compared them to see which assembly had the best output for the plasmid sequences of interest as described on my thesis [section 2.2.1].

#### Tool 01: Unicycler hybrid assembly
#### Input Data

| File | Format | Description |
|------|--------|-------------|
| input_file.fastq.gz | fastq.gz | [Nanopore long-reads] |
| input_file_R1.fastq.gz | fastq.gz | [Illumina short-reads R1] |
| input_file_R2.fastq.gz | fastq.gz | [Illumina short-reads R2] |

After obtaining all Illumina short-reads from the NCBI Genomes database and have them on the same folder as the filtered Nanopore long-reads, I prepared a text file and saved as "samples_list_unicycler.txt". This list is a tab separated file with three columns. The first column needs to be the name of the nanopore reads file, the second column the corresponding illumina pair 1 and the third column the corresponding Illumina pair 2. Example:

| Nanopore reads | Illumina reads R1 | Illumina reads R2 |
|------|--------|-------------|
| SC0119_FL.fastq.gz | SRR17899236_1.fastq.gz | SRR17899236_2.fastq.gz |
| SC0134_FL.fastq.gz | SRR17899220_1.fastq.gz | SRR17899220_2.fastq.gz |
| SC0212_FL.fastq.gz | SRR17900509_1.fastq.gz | SRR17900509_2.fastq.gz |

Then, run the Unicycler hybrid assembly script:

### Unicycler hybrid assembly script

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=hybrid_assembly_campy # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=168:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=4 # number of cpus/threads #usually I put 6cpus
#SBATCH --mem=20G            # RAM 

#access working directory
cd /nesi/nobackup/massey03742/Campy_SACNZ/reads/filtered_reads

module purge
module load Unicycler/0.5.0-gimkl-2020a-Python-3.8.2

# Path to the file containing the list of samples
SAMPLE_LIST="sample_list_unicycler.txt"

# Loop through each line in the sample list file
while IFS=$'\t' read -r nanopore_file illumina_file_1 illumina_file_2; do
    # Extract the base name for the output directory
    base_name=$(basename "${nanopore_file}" ".fastq.gz")

    # Create an output directory for this sample
    output_dir="/nesi/nobackup/massey03742/Campy_SACNZ/${base_name}"

    # Run Unicycler for the current sample
    unicycler -1 "${illumina_file_1}" -2 "${illumina_file_2}" -l "${nanopore_file}" -o "${output_dir}"

    echo "Finished assembly for ${base_name}"

done < "${SAMPLE_LIST}"

```
After the run is over, I extract the fasta files from each subdirectory that Unicycler creates for each genome and rename it with the proper name of the isolate, which is the name of the parent folder. I do that because inside each folder created by Unicycler, the name of the assembly is usually "assembly.fasta". So renaming the fasta file before moving to a directory where polishing will happen is important.

### Moving and renaming files

```bash

#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=moving_n_renaming_campy # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=168:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=2
#SBATCH --mem=500mb

cd /nesi/nobackup/massey03742/Campy_SACNZ/reads/filtered_reads

# Set the main directory
main_directory="/nesi/nobackup/massey03742/Campy_SACNZ/assemblies/"

# Loop through each directory named "SCXXXX"
for dir in "$main_directory"/SC*; do
    if [ -d "$dir" ]; then
        # Get the directory name (basename)
        dir_name=$(basename "$dir")
        
        # Define the old and new file names
        old_file="$dir/assembly.fasta"
        new_file="$main_directory/${dir_name}_unicycler.fasta"
        
        # Check if the assembled.fasta file exists
        if [ -f "$old_file" ]; then
            # Copy and rename the file to the main directory
            cp "$old_file" "$new_file"
            echo "Copied and renamed $old_file to $new_file"
        else
            echo "No assembled.fasta file found in $dir"
        fi
    fi
done  
```

#### Tool 02: Flye assembly
#### Input Data

| File | Format | Description |
|------|--------|-------------|
| input_file.fastq.gz | fastq.gz | [Nanopore long-reads] |


I also used Flye as an assembly tool, in order to compare both assemblies. So, before start polishing, I will show the script used for Flye. Then all resulting assemblies will follow the same polishing steps.

#### Flye assembly script

```bash

#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=flye2_assembly_campy # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=168:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=6 # number of cpus/threads, stick to 2 for now
#SBATCH --mem=15G            # RAM

cd /nesi/nobackup/massey03742/Campy_SACNZ/reads/filtered_reads

module purge
module load Flye/2.9.5-foss-2023a-Python-3.11.6

# Loop through each QC fastq.gz file in the current directory
for input_file in *.fastq.gz; do
    # Get the base name of the file (without the _FL.fastq.gz extension)
    base_name=$(basename "$input_file" _FL.fastq.gz)
    # Define the output directory
    output_dir="${base_name}_flye"
    
    echo "Starting assembly for $input_file..."
    
    # Check if the output directory exists
    if [ -d "$output_dir" ]; then
        echo "Output directory $output_dir exists. Resuming assembly..."
        flye --nano-hq "$input_file" -t 6 --resume --out-dir "$output_dir"
    else
        echo "No existing output directory. Starting a new assembly for $input_file..."
        flye --nano-hq "$input_file" -t 6 --out-dir "$output_dir"
    fi
    
    echo "Assembly for $input_file completed successfully."
done

echo "All assemblies completed."
```

#### Moving and renaming files

```bash

#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=moving_n_renaming_campy # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=168:00:00 # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=2
#SBATCH --mem=500mb

# Set the main directory
main_directory="/nesi/nobackup/massey03742/Campy_SACNZ/assemblies"

# Loop through each directory named "SCXXXX"
for dir in "$main_directory"/SC*; do
    if [ -d "$dir" ]; then
        # Get the directory name (basename)
        dir_name=$(basename "$dir")
        
        # Define the old and new file names
        old_file="$dir/assembly.fasta"
        new_file="$main_directory/${dir_name}.fasta"
        
        # Check if the assembled.fasta file exists
        if [ -f "$old_file" ]; then
            # Copy and rename the file to the main directory
            cp "$old_file" "$new_file"
            echo "Copied and renamed $old_file to $new_file"
        else
            echo "No assembled.fasta file found in $dir"
        fi
    fi
done   

```

## ENDS HERE ##

## Workflow Summary

[1-2 paragraph description of the analysis approach, key steps, and outputs]


### Input Data

| File | Format | Description |
|------|--------|-------------|
| input_file.fasta | FASTA | [Description] |
| metadata.csv | CSV | [Description] |
| reference.db | Database | [Description] |



## Usage

### Quick Start

```bash
# Navigate to this directory
cd [subsection-path]

# Review the script
cat script_name.sh

# Run the workflow (local)
bash script_name.sh

# OR submit to HPC
sbatch script_name.sbatch
```

### Step-by-Step Execution

#### Step 1: [First step name]

[Description of what happens and why]

```bash
bash step1_script.sh -i input.fasta -o output_dir
```

**Expected output:**
- `file_description.txt`
- `results_summary.csv`

#### Step 2: [Second step name]

[Description]

```bash
bash step2_script.sh -i output_dir/previous_output.txt
```

**Expected output:**
- `processed_results.xlsx`

#### Step 3: [Third step name]

[Description]

```bash
Rscript step3_analysis.R -i processed_results.xlsx -o figures/
```

**Expected output:**
- `figure1.pdf`
- `summary_stats.txt`

## Parameter Explanations

### Script: [script_name.sh]

```bash
-i, --input       Input FASTA file [required]
-o, --output      Output directory [required, default: ./output]
-t, --threads     Number of CPU threads [default: 4]
-m, --memory      Memory allocation in GB [default: 16]
-p, --prefix      Prefix for output files [default: analysis]
-v, --verbose     Enable verbose output [flag]
```

### Example command:
```bash
bash script_name.sh \
  -i input_sequences.fasta \
  -o ./results \
  -t 8 \
  -m 32 \
  -p my_analysis \
  -v
```

## Outputs

### Main Results

| File | Format | Contents |
|------|--------|----------|
| main_output.csv | CSV | [Description of contents] |
| analysis_results.xlsx | Excel | [Description] |
| summary_report.txt | Text | [Description] |

### Expected Output Summary

- Number of sequences processed: [X]
- Number of hits/matches: [Y]
- Main output file size: ~[X] MB
- Expected runtime: [X-Y] hours

## References and Documentation

- [Software name] documentation: [link]
- [Algorithm/method] paper: [citation]
- Related analysis: [Link to related subsection in thesis]

## Author Notes

[Any specific notes about implementation, deviations from standard methods, or important considerations]

---

*Last updated: [Date]*
*Tested on: [OS/Platform] with [Software versions]*
