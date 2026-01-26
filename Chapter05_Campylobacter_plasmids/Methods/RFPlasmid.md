# RFPlasmid Analysis

## Overview
Brief explanation of what RFPlasmid does and why you used it.

## Input Data
- 245 Campylobacter isolates from NZ (2005-2024)
- Assembled genomes in FASTA format
- Located at: `/path/to/assemblies/`

## Code
```bash
# Activate conda environment
conda activate rfplasmid

# Run RFPlasmid on all assemblies
for file in assemblies/*.fasta; do
    rfplasmid --input "$file" --output results/
done
```

## Output Interpretation
Explain what the output columns mean, how you filtered results, etc.

## Notes
Any issues you encountered, parameter choices, etc.
