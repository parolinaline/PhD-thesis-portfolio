# [Subsection Title] Methodology

## Overview

[Brief description of what this subsection does and its purpose in the larger analysis workflow]

## Workflow Summary

[1-2 paragraph description of the analysis approach, key steps, and outputs]

## Requirements

### Software and Packages

- [Software name] (v.x.x) - [purpose]
- [Software name] (v.x.x) - [purpose]

**R packages:**
```R
# Install if needed
install.packages(c("package1", "package2"))
# Or from Bioconductor:
BiocManager::install(c("package1"))
```

**Python packages:**
```bash
pip install package1 package2
```

### Input Data

| File | Format | Description |
|------|--------|-------------|
| input_file.fasta | FASTA | [Description] |
| metadata.csv | CSV | [Description] |
| reference.db | Database | [Description] |

### Computing Resources

- **CPU cores**: X (per job)
- **Memory**: X GB
- **Time estimate**: X hours
- **HPC platform**: NeSI / Local
- **Storage**: X GB (temporary)

## Files in This Directory

```
├── script_name.sh/.py/.R    # Main workflow script
├── script_helper.py         # Helper functions (if applicable)
├── README.md                # This file
├── input/                   # Input data (or symlinks to data)
│   ├── raw/
│   └── processed/
└── output/                  # Analysis outputs
    ├── results/
    └── logs/
```

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

## Troubleshooting

### Common Issues

**Error: "Command not found: [software]"**
- Ensure software is installed: `which software_name`
- Check PATH: `echo $PATH`
- Load module on HPC: `module load software_name`

**Error: "Insufficient memory"**
- Increase memory allocation: `-m 64` (for 64 GB)
- Reduce number of threads: `-t 4`
- Process data in smaller batches

**Error: "Permission denied"**
- Check file permissions: `ls -la`
- Make script executable: `chmod +x script_name.sh`

**Output files not generated**
- Check error messages: `cat logs/error.log`
- Verify input file format: `head input_file.fasta`
- Run script with verbose flag: `-v`

### Getting Help

If you encounter issues:
1. Check the error log: `logs/error.log` or SLURM output file
2. Review input file format and contents
3. Test with a smaller subset of data
4. Review software documentation

## Outputs

### Main Results

| File | Format | Contents |
|------|--------|----------|
| main_output.csv | CSV | [Description of contents] |
| analysis_results.xlsx | Excel | [Description] |
| summary_report.txt | Text | [Description] |

### Supporting Files

- `logs/` - Workflow execution logs
- `temp/` - Intermediate files (can be deleted after verification)
- `figures/` - Generated plots and visualizations

## Verification

### Check Your Results

Run this verification script to ensure outputs are correct:

```bash
bash verify_output.sh output_directory/
```

Or manually verify:

```bash
# Check file sizes are reasonable
ls -lh output/*.csv output/*.xlsx

# Quick data inspection
head -20 output/main_output.csv
wc -l output/main_output.csv  # Count rows

# Check for errors in logs
grep -i "error\|warning" logs/*.log
```

### Expected Output Summary

- Number of sequences processed: [X]
- Number of hits/matches: [Y]
- Main output file size: ~[X] MB
- Expected runtime: [X-Y] hours

## Performance Notes

### Runtime Estimates

- Local execution (4 cores): ~X hours
- NeSI execution (8 cores): ~X hours
- Scales roughly: O(n log n) with data size

### Memory Usage

- Peak memory: ~X GB (for typical dataset)
- Memory scales linearly with data size

### Tips for Large Datasets

- Increase number of threads: `-t 16` or higher
- Use parallel processing: `GNU parallel` or `xargs`
- Process in chunks: see `split_data.sh`
- Check NeSI job parameters in `.sbatch` file

## References and Documentation

- [Software name] documentation: [link]
- [Algorithm/method] paper: [citation]
- Related analysis: [Link to related subsection in thesis]

## Author Notes

[Any specific notes about implementation, deviations from standard methods, or important considerations]

---

*Last updated: [Date]*
*Tested on: [OS/Platform] with [Software versions]*
