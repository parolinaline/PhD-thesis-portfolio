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

## Repository Structure

```
thesis-portfolio/
├── index.md                    # Home page
├── chapters.md                 # Chapter overview and navigation
├── about.md                    # Project background
├── _config.yml                 # Jekyll configuration
├── _layouts/
│   └── default.html           # Custom page template
├── assets/
│   └── css/
│       └── style.css          # Dark theme styling
└── chapters/
    ├── 01-introduction/
    ├── 02-literature/
    ├── 03-campylobacter-plasmids/
    │   ├── index.md
    │   ├── 3.1-background/
    │   ├── 3.2-methodology/
    │   │   ├── 3.2.1-assembly/
    │   │   ├── 3.2.2-plasmid-class/
    │   │   └── 3.2.3-phylogenetics/
    │   ├── 3.3-results/
    │   │   ├── 3.3.1-plasmid-diversity/
    │   │   ├── 3.3.2-comparative-genomics/
    │   │   └── 3.3.3-phylogenetics/
    │   └── 3.4-data/
    ├── 04-salmonella-surveillance/
    │   ├── index.md
    │   ├── 4.1-background/
    │   ├── 4.2-methodology/
    │   │   ├── 4.2.1-cgmlst-setup/
    │   │   ├── 4.2.2-mlst-analysis/
    │   │   ├── 4.2.3-amr-profiling/
    │   │   ├── 4.2.4-phylogenetics/
    │   │   └── 4.2.5-plasmids/
    │   ├── 4.3-results/
    │   │   ├── 4.3.1-sequence-types/
    │   │   ├── 4.3.2-amr-analysis/
    │   │   ├── 4.3.3-plasmids/
    │   │   └── 4.3.4-comparative-genomics/
    │   └── 4.4-data/
    ├── 05-comparative-analysis/
    └── 06-discussion/
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
- Data dictionary files explaining table contents
- Links between related sections

### Data

- Raw sequencing data (assemblies in FASTA format)
- Processed intermediate files
- Final results tables (Excel, CSV, TSV)
- Metadata files with sample information

## Using This Repository

### For Thesis Reviewers

1. **Start with [Chapters](chapters.html)** to navigate by topic
2. **Follow the methodology section** to understand analysis steps
3. **View results sections** for key findings
4. **Access raw/processed data** to verify outputs
5. **Review scripts** for implementation details

### For Reproducibility

Each workflow includes:
- Complete script code
- Input data files (or links to larger repositories)
- Step-by-step execution order
- Expected outputs for verification

To run analyses:
```bash
cd chapters/[chapter]/[subsection]/
cat README.md          # Review the method
bash script_name.sh    # Run the workflow
```

### For R Visualizations

All R scripts use common packages:
- `ggplot2` - Publication-quality graphics
- `ggtree` - Phylogenetic tree visualization
- `gggenomes` - Comparative genomics
- `dplyr`, `tidyr` - Data manipulation
- Standard statistics packages

Run R scripts:
```bash
Rscript script_name.R
```

## Theme and Design

This site uses:
- **Dark theme** for reduced eye strain during extended reading
- **Minimal design** with focus on content and readability
- **Clean typography** for academic clarity
- **Direct code viewing** in the browser
- **Responsive layout** for desktop and mobile viewing

Custom CSS styling in `assets/css/style.css`.

## Making Changes

### Adding a New Chapter

1. Create a new folder: `chapters/XX-chapter-name/`
2. Create `index.md` with content following existing templates
3. Create subfolders for methodology, results, data
4. Add links to `chapters.md`

### Adding Scripts

Place scripts in appropriate methodology folders:
```
chapters/XX/XX.Y-subsection/
├── script_name.sh/.py/.R
├── README.md
└── data/
    ├── input/
    └── output/
```

### Updating the Site

The site rebuilds automatically when you push to GitHub. To test locally:
```bash
bundle exec jekyll serve
```

## Privacy and Access Control

### Current Status
**Private repository** - Access restricted to specified reviewers

### Before Defense
- Keep repository private
- Share access links with supervisors/committee

### After Defense
- Make repository public
- Archive on Zenodo or institutional repository
- Update README with DOI and citation information

## Contact and Citation

For questions about this research, contact through your institution.

**Citation format** (to be updated after defense):
```
Author Name (2025). PhD Thesis Supplementary Materials. 
[Thesis title]. Repository: [GitHub URL]
```

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

## License

[Specify your license - e.g., CC BY 4.0 for academic attribution]

---

*Last updated: November 2025*

For more information, visit the [home page](index.md) or [chapters overview](chapters.md).
