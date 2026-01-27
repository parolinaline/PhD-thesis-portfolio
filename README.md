# PhD Thesis Supplementary Materials

A comprehensive, organized repository of scripts, workflows, data, and documentation supporting a PhD thesis on bacterial genomics and antimicrobial resistance in *Campylobacter* and *Salmonella*.

## Overview

This repository serves as supplementary materials for thesis chapters, enabling reviewers and readers to:
- Access complete analysis workflows and scripts
- Review intermediate data files and processing steps
- Examine results tables and supplementary figures
- Understand methodology in detail
- Reproduce analyses

## Quick Start

### View the Site

This repository is published as a GitHub Pages website. Visit:
```
https://username.github.io/thesis-portfolio/
```

### Local Setup

If you want to build and test locally:

```bash
# Install Jekyll (if needed)
gem install bundler jekyll

# Clone the repository
git clone https://github.com/username/thesis-portfolio.git
cd thesis-portfolio

# Build and serve locally
bundle exec jekyll serve

# Visit http://localhost:4000/
```


## Chapter Organization

Each chapter follows a consistent structure:

- **Introduction/Background**: Context and research questions
- **Methodology**: Step-by-step workflows with documented scripts
  - Bash scripts for HPC pipelines
  - Python scripts for complex analysis
  - R scripts for visualization and statistics
- **Results**: Analysis outputs, figures, and supplementary tables
- **Data**: Raw and processed data files

## Key Features

### Scripts

All scripts are documented with:
- Clear comments explaining each step
- Usage examples and parameter descriptions
- Dependencies and software versions
- Expected input/output formats

**Languages**: R, Bash, Python

### Workflows

HPC-optimized workflows using:
- NeSI SLURM/SBATCH job scheduling
- Parallel processing for large-scale analyses
- Standard bioinformatics tools (chewBACCA, GrapeTree, BLAST, etc.)

### Documentation

Each section includes:
- README files with detailed methods
- Step-by-step guides for complex analyses
- Links between related sections

### Data

- Raw sequencing data (assemblies in FASTA format)
- Processed intermediate files
- Final results tables (Excel, CSV, TSV)
- Metadata files with sample information

## Using This Repository

### For Thesis Reviewers

1. **Start with [Chapters](./chapters.md)** to navigate by topic
2. **Follow the methodology section** to understand analysis steps
3. **View results sections** for key findings
4. **Access raw/processed data** to verify outputs
5. **Review scripts** for implementation details


### For R Visualizations

All R scripts use common packages:
- `ggplot2` - Publication-quality graphics
- `ggtree` - Phylogenetic tree visualization
- `gggenomes` - Comparative genomics
- `dplyr`, `tidyr` - Data manipulation
- Standard statistics packages


## Privacy and Access Control

### Current Status
**Private repository** - Access restricted to specified reviewers


---

## Tools and Software Used

### Bioinformatics
- **PubMLST** - MLST/cgMLST schemes
- **chewBACCA** - cgMLST allele calling
- **GrapeTree** - Phylogenetic visualization
- **BLAST** - Sequence searching
- **MOB-suite** - Plasmid analysis

### Programming Languages
- **R** (v4.x) - Phylogenetics, visualization, statistics
- **Python** (v3.x) - Bioinformatics pipelines, data processing
- **Bash** - HPC workflow automation

### R Packages
```R
ggtree, aplot, gggenomes, ggplot2, dplyr, tidyr
```

### Computing Infrastructure
- **NeSI** (New Zealand eScience Infrastructure)
- SLURM job scheduler
- High-performance computing for large-scale analyses


---

*Last updated: November 2025*
