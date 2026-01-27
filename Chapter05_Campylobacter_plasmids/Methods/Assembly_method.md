# Assemblying genomes

## Overview

Here I explain how the assemblies were obtained for Chapter 05, starting from the filtering process of the reads until the polishing steps. I used two different approaches, and the final assemblies were compared and chosen as described on the thesis. The tools and versions used are indicated in the scripts. 

## Workflow Summary

1 - Filtering Nanopore long-reads using FiltLong
2 - Assemblying genomes:
    a- Unicycler hybrid assembly mode using Nanopore long-reads and Illumina short-reads
    b- Flye assembly using only Nanopore long-reads
3- Polish assemblies using Nanopore long-reads and Illumina short-reads 


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
cd /nesi/nobackup/massey03742/Auckland_results/super_high_accuracy_fastq_files/campy

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


COPY FROM HERE
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
