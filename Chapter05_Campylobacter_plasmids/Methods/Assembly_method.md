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
### Step 1: [Filtering Nanopore long-reads after basecalling]
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

### Step 2: [Assemblying genomes]

I had the Illumina short-reads available for these isolates: BioProject PRJNA675916
I tested two assembly tools and compared them to see which assembly had the best output for the plasmid sequences of interest as described on my thesis [section 2.2.1].

### Tool 01: Unicycler hybrid assembly
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

#### Unicycler hybrid assembly script

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

#### Moving and renaming files

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

### Tool 02: Flye assembly
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
### Step 3: [Polishing with Racon]

After obtaining all assemblies from both assembly tools (Flye and Racon), let's start the polishing proccess with Racon. Racon uses long-reads for that. Before starting to run Racon, let's prepare two text files. One with the name and path of the reads that will be used, and the other one with the name of the assemblies that will be used. The read that correspond to a certain assembly need to be on the same row number of the assembly. For example:
File "long_reads.txt"
| File | 
|------|
| path/to/filtered/reads/SC0119_FL.fastq.gz | 
| path/to/filtered/reads/SC0134_FL.fastq.gz | 
| path/to/filtered/reads/SC0212_FL.fastq.gz |
| path/to/filtered/reads/SC0119_FL.fastq.gz | 
| path/to/filtered/reads/SC0134_FL.fastq.gz | 
| path/to/filtered/reads/SC0212_FL.fastq.gz |


File "Initial_assemblies.txt"
| File | 
|------|
| path/to/assembly/SC0119_unicycler.fasta | 
| path/to/assembly/SC0134_unicycler.fasta | 
| path/to/assembly/SC0212_unicycler.fasta |
| path/to/assembly/SC0119_flye.fasta | 
| path/to/assembly/SC0134_flye.fasta | 
| path/to/assembly/SC0212_flye.fasta |

The script will connect the assembly to the reads based on the row of the "long_reads.txt" and "Initial_assemblies.txt". I ran 4 iterations of the script following the methodology described here:
<https://doi.org/10.1186/s12864-024-10582-x>

#### Input Data

| File | Format | Description |
|------|--------|-------------|
| input_file.fastq.gz | fastq.gz | [Nanopore long-read already filtered] |
| input_file.fasta | fasta | Assembly |

#### Racon script

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=campy_racon_massey03742
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00


module purge 
module load Racon/1.5.0-GCC-11.3.0
module load minimap2/2.28-GCC-12.3.0

# Parameters
threads=4
maxiter=4
long_reads_file="long_reads.txt"
initial_assemblies_file="initial_assemblies.txt"

# Validate input files
if [[ ! -f ${long_reads_file} ]]; then
    echo "Error: '${long_reads_file}' not found. Please provide a valid file."
    exit 1
fi

if [[ ! -f ${initial_assemblies_file} ]]; then
    echo "Error: '${initial_assemblies_file}' not found. Please provide a valid file."
    exit 1
fi

# Process both files line by line
paste "${long_reads_file}" "${initial_assemblies_file}" | while IFS=$'\t' read -r long_reads initial_assembly; do
    isolate_id=$(basename ${long_reads} _FL.fastq.gz)
    output_dir="/nesi/nobackup/massey03742/Campy_SACNZ/racon/${isolate_id}"

    # Validate files exist
    if [[ ! -f ${long_reads} ]]; then
        echo "Error: Long read file '${long_reads}' not found. Skipping..."
        continue
    fi
    if [[ ! -f ${initial_assembly} ]]; then
        echo "Error: Initial assembly file '${initial_assembly}' not found. Skipping..."
        continue
    fi

    # Create output directory
    mkdir -p ${output_dir} || { echo "Failed to create output directory ${output_dir}"; exit 1; }

    # Perform Racon polishing for maxiter iterations
    polished=${initial_assembly}
    for iter in $(seq 1 ${maxiter}); do
        minimap2 -x map-ont -t ${threads} ${polished} ${long_reads} > ${output_dir}/mappings_${isolate_id}_${iter}.paf || { echo "Error during minimap2. Skipping..."; break; }
        racon -t ${threads} ${long_reads} ${output_dir}/mappings_${isolate_id}_${iter}.paf ${polished} > ${output_dir}/racon_round_${isolate_id}_${iter}.fasta || { echo "Error during racon. Skipping..."; break; }
        polished=${output_dir}/racon_round_${isolate_id}_${iter}.fasta
    done

    echo "Finished Racon for ${isolate_id}"
done
```
Similarly to what was done with Unicycler and Flye, I extracted the last fasta file (round 4) from Racon results from the subdirectories of each isolate and move them to one specific directory in order to run Medaka.
#### Moving Racon files
```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=copying_files_racon_massey03742
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00


# Set the parent directory containing the isolate directories
parent_dir="/nesi/nobackup/Campy_SACNZ/Racon"

# Set the output directory where the final Racon-polished files will be copied
output_dir="/nesi/nobackup/Campy_SACNZ/medaka"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop through each isolate directory
for isolate_dir in "$parent_dir"/*; do
    if [[ -d "$isolate_dir" ]]; then
        # Get the isolate ID from the directory name
        isolate_id=$(basename "$isolate_dir")

        # Find the latest Racon-polished file in the isolate directory
        latest_file=$(find "$isolate_dir" -type f -name "racon_round_${isolate_id}_*.fasta" | sort | tail -n 1)

        if [[ -n "$latest_file" ]]; then
            # Copy the file to the output directory
            cp "$latest_file" "$output_dir/"

            echo "Copied: $latest_file to $output_dir"
        else
            echo "No Racon-polished file found for isolate: $isolate_id"
        fi
    fi
done

echo "All files copied to $output_dir."
```

### Step 4: [Polishing with Medaka]

After obtaining all assemblies from Racon, I procceded with the second polishing step using Medaka. The input file for medaka is a file called "sample_files.txt". This is a tab delimited file. This file needs to be organized as following:
| Sample_ID | Path_to_Racon_file| Path_to_long_reads  
|------|------|------|
| SC0119 | path/to/racon/SC0119_4.fasta | path/to/filtered/reads/SC0119_FL.fastq.gz |

#### Medaka script

```bash

#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=medaka_campy_massey03742 # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=15G            # RAM
#SBATCH --cpus-per-task=4    # number of cpus/threads
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)



# Load necessary modules
module purge
module load medaka/1.11.1-Miniconda3-22.11.1-1
module load BWA/0.7.17-gimkl-2017a
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0

# Set parameters
input_file="sample_files.txt" # Path to the input text file
output_dir="/nesi/nobackup/massey03742/Campy_SACNZ/medaka"
threads=4

# Ensure the input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

# Read the text file line by line
while IFS=$'\t' read -r sample_id fasta_file long_read; do
    if [[ -z "$sample_id" || -z "$fasta_file" || -z "$long_read" ]]; then
        echo "Error: Missing data in line: $sample_id, $fasta_file, $long_read"
        continue
    fi

    # Check if the files exist
    if [[ ! -f "$fasta_file" ]]; then
        echo "Error: FASTA file '$fasta_file' not found for sample '$sample_id'!"
        continue
    fi
    if [[ ! -f "$long_read" ]]; then
        echo "Error: Long read file '$long_read' not found for sample '$sample_id'!"
        continue
    fi

    # Index the FASTA file
    echo "Indexing FASTA file: $fasta_file"
    samtools faidx "$fasta_file"

    # Run Medaka consensus
    echo "Running Medaka for sample: $sample_id"
    medaka_consensus -i "$long_read" -d "$fasta_file" -o "${output_dir}/${sample_id}_medaka" -t "$threads"

    # Move final polished assembly
    if [[ -f "${output_dir}/${sample_id}_medaka/consensus.fasta" ]]; then
        mv "${output_dir}/${sample_id}_medaka/consensus.fasta" "${output_dir}/${sample_id}_medaka.fasta"
    else
        echo "Warning: Medaka consensus failed for sample '$sample_id'. No consensus.fasta found."
    fi
done < "$input_file"

```
### Step 5: [Polishing with NextPolish]

After obtaining all polished assemblies from Medaka, I procceded with the third polishing step using NextPolish. NextPolish uses short-reads to polish the assemblies. Before running the NextPolish script, I prepared the alignment files that are used in NextPolish. First, I prepared a file called "files_list.txt" with the path for the medaka fasta file and the corresponding paired short-reads.   This is a tab delimited file. This file needs to be organized as following:

| Sample_ID | Path_to_Medaka_file| Path_to_short_reads_R1 | Path_to_short_reads_R2 
|------|------|------|------|
| SC0119 | path/to/racon/SC0119_medaka.fasta | path/to/filtered/reads/SC0119_R1.fastq.gz | path/to/filtered/reads/SC0119_R2.fastq.gz |

#### Preparing isolates before running NextPolish

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=preparetonextpolish_massey03742 # job name
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=10G            # RAM
#SBATCH --cpus-per-task=8     # number of cpus/threads
#SBATCH --time=10:00:00       # Walltime (HH:MM:SS)


module purge
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0

# Parameters
threads=8
input_dir="/nesi/nobackup/massey03742/Campy_SACNZ/medaka/"
output_dir="/nesi/nobackup/massey03742/Campy_SACNZ/medaka/nextpolish"

# Ensure the output directory exists
mkdir -p ${output_dir}

# Check if the file list exists
file_list="files_list.txt"
if [[ ! -f ${file_list} ]]; then
    echo "Error: File list '${file_list}' not found. Please provide a valid file list."
    exit 1
fi

# Loop through each line in the file list; it must be tab-separated
while read -r isolate fasta read1 read2; do
    echo "Processing $isolate..."

    # Validate input files
    if [[ ! -f ${input_dir}${fasta} ]]; then
        echo "Error: FASTA file '${input_dir}${fasta}' not found for isolate $isolate. Skipping..."
        continue
    fi
    if [[ ! -f ${input_dir}${read1} || ! -f ${input_dir}${read2} ]]; then
        echo "Error: Read files '${input_dir}${read1}' and/or '${input_dir}${read2}' not found for isolate $isolate. Skipping..."
        continue
    fi

    # Ensure the reference genome is indexed
    samtools faidx ${input_dir}${fasta}
    echo "Index of $fasta done"

    # Align reads with minimap2
    minimap2 -ax sr -t ${threads} ${input_dir}${fasta} ${input_dir}${read1} ${input_dir}${read2} > ${output_dir}${isolate}.sam || { echo "Minimap2 alignment failed for $isolate. Skipping..."; continue; }
    echo "Alignment for $isolate done"

    # Convert SAM to BAM, filter out unmapped reads, and fix mate information
    samtools view -@ ${threads} -F 0x4 -b ${output_dir}${isolate}.sam | \
        samtools fixmate -m --threads ${threads} - - | \
        samtools sort -@ ${threads} -o ${output_dir}${isolate}.sorted.bam || { echo "SAM to BAM conversion failed for $isolate. Skipping..."; continue; }
    echo "BAM file for $isolate ready"

    # Mark duplicates
    samtools markdup --threads ${threads} -r ${output_dir}${isolate}.sorted.bam ${output_dir}${isolate}.sorted.markdup.bam || { echo "Duplicate marking failed for $isolate. Skipping..."; continue; }
    echo "Duplicate marking for $isolate done"

    # Index the marked duplicate BAM file
    samtools index -@ ${threads} ${output_dir}${isolate}.sorted.markdup.bam || { echo "Indexing failed for $isolate. Skipping..."; continue; }
    echo "Indexing for $isolate done"

    # Cleanup SAM file to save space
    rm -f ${output_dir}${isolate}.sam
    echo "Temporary SAM file for $isolate removed"

done < ${file_list}

echo "All isolates processed successfully."

```

#### NextPolish

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=preparetonextpolish_massey03742 # job name
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=10G            # RAM
#SBATCH --cpus-per-task=8     # number of cpus/threads
#SBATCH --time=10:00:00       # Walltime (HH:MM:SS)


module purge
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0

# Parameters
threads=8
input_dir="/nesi/nobackup/massey03742/Campy_SACNZ/medaka/nextpolish"
output_dir="/nesi/nobackup/massey03742/Campy_SACNZ/medaka/nextpolish"

# Ensure the output directory exists
mkdir -p ${output_dir}

# Check if the file list exists
file_list="files_list.txt"
if [[ ! -f ${file_list} ]]; then
    echo "Error: File list '${file_list}' not found. Please provide a valid file list."
    exit 1
fi

# Loop through each line in the file list; it must be tab-separated
while read -r isolate fasta read1 read2; do
    echo "Processing $isolate..."

    # Validate input files
    if [[ ! -f ${input_dir}${fasta} ]]; then
        echo "Error: FASTA file '${input_dir}${fasta}' not found for isolate $isolate. Skipping..."
        continue
    fi
    if [[ ! -f ${input_dir}${read1} || ! -f ${input_dir}${read2} ]]; then
        echo "Error: Read files '${input_dir}${read1}' and/or '${input_dir}${read2}' not found for isolate $isolate. Skipping..."
        continue
    fi

    # Ensure the reference genome is indexed
    samtools faidx ${input_dir}${fasta}
    echo "Index of $fasta done"

    # Align reads with minimap2
    minimap2 -ax sr -t ${threads} ${input_dir}${fasta} ${input_dir}${read1} ${input_dir}${read2} > ${output_dir}${isolate}.sam || { echo "Minimap2 alignment failed for $isolate. Skipping..."; continue; }
    echo "Alignment for $isolate done"

    # Convert SAM to BAM, filter out unmapped reads, and fix mate information
    samtools view -@ ${threads} -F 0x4 -b ${output_dir}${isolate}.sam | \
        samtools fixmate -m --threads ${threads} - - | \
        samtools sort -@ ${threads} -o ${output_dir}${isolate}.sorted.bam || { echo "SAM to BAM conversion failed for $isolate. Skipping..."; continue; }
    echo "BAM file for $isolate ready"

    # Mark duplicates
    samtools markdup --threads ${threads} -r ${output_dir}${isolate}.sorted.bam ${output_dir}${isolate}.sorted.markdup.bam || { echo "Duplicate marking failed for $isolate. Skipping..."; continue; }
    echo "Duplicate marking for $isolate done"

    # Index the marked duplicate BAM file
    samtools index -@ ${threads} ${output_dir}${isolate}.sorted.markdup.bam || { echo "Indexing failed for $isolate. Skipping..."; continue; }
    echo "Indexing for $isolate done"

    # Cleanup SAM file to save space
    rm -f ${output_dir}${isolate}.sam
    echo "Temporary SAM file for $isolate removed"

done < ${file_list}

echo "All isolates processed successfully."


```

All fasta files can now be moved to the desired directory to start in-silico analyses.

## ENDS HERE ##


---

*Last updated: [Date]*
*Tested on: [OS/Platform] with [Software versions]*
