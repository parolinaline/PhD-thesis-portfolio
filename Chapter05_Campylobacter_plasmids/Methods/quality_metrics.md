# Quality metrics of Nanopore reads and hybrid assemblies

## Overview
Here I explain all tools used to assess the quality of reads and hybrid assemblies used for Chapter 05.

## Workflow Summary
1. Using Nanoplot to assess read quality
2. Tidy up the data and generate a summary file
3. Rename NanoStat files before running MultiQC
4. MultiQC to visualise a summary of reads quality
5. Run QUAST to assess assemblies quality

### Step-by-Step Execution
### Step 1: [Running Nanoplot]
#### Script
```bash
#!/bin/bash -e
#SBATCH --job-name=nanoplot_campy
#SBATCH --output=slurmlog/%x.%j_16jan24.out
#SBATCH --error=slurmlog/%x.%j_16jan24.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

module purge
module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6

# Directory with your fastq.gz files
READS_DIR="/nesi/nobackup/massey03742/Campy_SACNZ/filtered_reads"
OUTPUT_DIR="/nesi/nobackup/massey03742/Campy_SACNZ/filtered_reads/nanoplot"

# Create main output folder
mkdir -p "$OUTPUT_DIR"

# Loop over all FASTQ.gz files
for file in "$READS_DIR"/*.fastq.gz; do
    # Extract base filename without extension
    base=$(basename "$file" .fastq.gz)

    # Make output subfolder for each sample
    outdir="$OUTPUT_DIR/$base"
    mkdir -p "$outdir"

    # Run NanoPlot
    NanoPlot --fastq "$file" -o "$outdir" --threads 4 --N50 --plots hex dot

    echo "NanoPlot done for $base"
done

```

### Step 2: [Generate a summary file from Nanoplot output]
### Overview
This script will extract the main information from the Nanoplot output and generate a summary file.
#### Script

```bash
#!/bin/bash -e
#SBATCH --job-name=nanoplot_summary_results
#SBATCH --output=slurmlog/%x.%j_16jan24.out
#SBATCH --error=slurmlog/%x.%j_16jan24.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

# Directory containing NanoPlot output folders
INPUT_DIR="/nesi/nobackup/massey03742/Campy_SACNZ/filtered_reads/nanoplot"
OUTPUT_TSV="nanoplot_summary.tsv"

# Write header
echo -e "ID\tmean_len\tmean_qual\tmedian_len\tmedian_qual\tn_reads\tN50\tstdev_len\ttotal_bases\tpct_Q10\tpct_Q15\tpct_Q20\tpct_Q25\tpct_Q30" > "$OUTPUT_TSV"

# Loop through each folder
for folder in "$INPUT_DIR"/*/; do
    stats_file="${folder}/NanoStats.txt"
    [ -f "$stats_file" ] || continue

    # Extract ID from folder name
    ID=$(basename "$folder")

    # Use awk or grep+cut to extract fields
    mean_len=$(grep "Mean read length:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    mean_qual=$(grep "Mean read quality:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    median_len=$(grep "Median read length:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    median_qual=$(grep "Median read quality:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    n_reads=$(grep "Number of reads:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    N50=$(grep "Read length N50:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    stdev_len=$(grep "STDEV read length:" "$stats_file" | awk -F ':' '{print $2}' | xargs)
    total_bases=$(grep "Total bases:" "$stats_file" | awk -F ':' '{print $2}' | xargs)

    pct_Q10=$(grep "^>Q10:" "$stats_file" | grep -oP '\(\K[0-9.]+(?=%\))')
    pct_Q15=$(grep "^>Q15:" "$stats_file" | grep -oP '\(\K[0-9.]+(?=%\))')
    pct_Q20=$(grep "^>Q20:" "$stats_file" | grep -oP '\(\K[0-9.]+(?=%\))')
    pct_Q25=$(grep "^>Q25:" "$stats_file" | grep -oP '\(\K[0-9.]+(?=%\))')
    pct_Q30=$(grep "^>Q30:" "$stats_file" | grep -oP '\(\K[0-9.]+(?=%\))')

    # Append to output TSV
    echo -e "${ID}\t${mean_len}\t${mean_qual}\t${median_len}\t${median_qual}\t${n_reads}\t${N50}\t${stdev_len}\t${total_bases}\t${pct_Q10}\t${pct_Q15}\t${pct_Q20}\t${pct_Q25}\t${pct_Q30}" >> "$OUTPUT_TSV"
done

echo "NanoPlot summary saved to $OUTPUT_TSV"
```

### Step 3: [Rename NanoStat files and run MultiQC]
### Overview
Nanoplot output contains one directory for each read used as input. Inside each folder you have all results from Nanoplot.
This script will rename all NanoStats.txt files inside each subdirectory according to the name of the isolate's ID, which in my case if the name of the parent directory,
and I moved all renamed NanoStats files to a separate directory to run MultiQC. This step is necessary otherwise MultiQC interprets all NanoStats.txt inside each isolate's folder from Nanoplot output
as the same file and overwrites it.

#### Script

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=collect_nanostats
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G

# Source directory with nanoplot subfolders
SOURCE_DIR="/nesi/nobackup/massey03742/Campy_SACNZ/filtered_reads/nanoplot"

# Destination directory for renamed NanoStats files
DEST_DIR="/nesi/nobackup/massey03742/Campy_SACNZ/filtered_reads/nanoplot/nanostats_collected"

# Create destination if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all NanoStats.txt files in subdirectories
for file in "$SOURCE_DIR"/*/NanoStats.txt; do
    # Skip if no matches found
    [[ -e "$file" ]] || continue
    
    # Get the parent directory name (isolate name) and path
    parent_dir=$(dirname "$file")
    isolate_name=$(basename "$parent_dir")
    
    # New filename: isolate_name_NanoStats.txt
    new_name="${isolate_name}_NanoStats.txt"
    
    # Rename in original location
    mv "$file" "$parent_dir/$new_name"
    echo "Renamed: $file -> $parent_dir/$new_name"
    
    # Copy to destination folder
    cp "$parent_dir/$new_name" "$DEST_DIR/$new_name"
    echo "Copied to: $DEST_DIR/$new_name"
done

echo "Done! Files renamed in place and copied to $DEST_DIR"
```
### Step 4: [Run MultiQC]
### Overview
MultiQC gets all NanoStats output from Nanoplot and summarises the main information with visual and clear information.

#### Script

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=multiqc # job name (shows up in the queue)
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G            # RAM
#SBATCH --cpus-per-task=2 # number of cpus/threads
#SBATCH --time=24:00:00 # Walltime (HH:MM:SS)

cd /nesi/nobackup/massey03742/Campy_SACNZ/filtered_reads/nanoplot/nanostats_collected

module purge
module load MultiQC/1.24.1-foss-2023a-Python-3.11.6


multiqc . -o multiqc_report_campy
```

### Step 5: [Run QUAST]
### Overview
Quast provides quality metrics for genome assemblies. Here I am using as input all polished assemblies generated with Nanopore long reads and Illumina short reads following what has been described in (./Assembly_method.md).

#### Script

```bash
#!/bin/bash -e
#SBATCH --job-name=quast_campy
#SBATCH --output=slurmlog/%x.%j_16jan24.out
#SBATCH --error=slurmlog/%x.%j_16jan24.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4 
#SBATCH --account=massey03742
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS)

module purge
module load QUAST/5.2.0-gimkl-2022a

# Directory containing your assemblies (FASTA files)
ASSEMBLY_DIR="/nesi/nobackup/massey03742/Campy_SACNZ/best_assemblies"
OUTPUT_DIR="quast_results"

# Create main output folder
mkdir -p "$OUTPUT_DIR"

# Loop through all .fasta or .fa files
for file in "$ASSEMBLY_DIR"/*.fasta "$ASSEMBLY_DIR"/*.fa; do
    # Skip if no files found
    [ -e "$file" ] || continue

    # Get base filename without extension
    base=$(basename "$file" .fasta)
    base=$(basename "$base" .fa)

    # Set output directory for this sample
    outdir="$OUTPUT_DIR/$base"
    mkdir -p "$outdir"

    # Run QUAST without reference
    quast.py "$file" -o "$outdir" --threads 4 --gene-finding

    echo "QUAST finished for $base"
done
```

The output files can be seen in Appendix (insert here the links to files).
