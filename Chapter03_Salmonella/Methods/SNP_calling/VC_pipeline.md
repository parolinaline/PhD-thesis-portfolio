# Methodology
# Calling SNPs using minimap2 and Clair3

## Overview

Here I explain the steps involved in generating a full alignment and a core SNP alignment using Nanopore long reads.
The mapping was obtained using minimap2 and the VCF files were generated using Clair3.

This pipeline replicates the Snippy workflow for long reads by:
1. Mapping reads to reference
2. Calling variants with Clair3 (haploid mode)
3. Creating callable position masks (positions with sufficient coverage)
4. Generating per-isolate consensus sequences with low-coverage positions masked as N
5. Using snp-sites to extract core alignments

## Workflow Summary

1. Mapping reads against a reference using minimap2
2. Call SNPs (generate VCF files) using Clair3
3. Filter VCF files using BCFtools
4. Generate callable position masks per isolate
5. Generate consensus FASTA per isolate
6. Extract core genome alignment using snp-sites
7. Filter recombination with Gubbins
8. Generate SNP distance matrix

---

## Step 1: Mapping reads against a reference

### Input Data

Prepare a text file with the name of all reads that will be used in this analysis. I am using Nanopore long-reads.
The file will be named `list_all.txt` and should look like this:

| Filename |
|----------|
| SGV002_FL.fastq.gz |
| SGV003_FL.fastq.gz |
| SGV004_FL.fastq.gz |
| SGV005_FL.fastq.gz |

Then we can run the script using minimap2 as the tool for mapping. This will generate BAM files necessary to run Clair3.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=minimap2_bovis
#SBATCH --output=slurmlog/%x.%A_%a.out
#SBATCH --error=slurmlog/%x.%A_%a.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --array=1-366%50

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/reads

module load SAMtools/1.22-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0

REF=SBV005_medaka.fasta

# Index reference only once (check if exists)
[[ ! -f "${REF}.fai" ]] && samtools faidx ${REF}

# Get the read file for this array task
read=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_all.txt)
base=$(basename "$read" .fastq.gz)

echo "Task ${SLURM_ARRAY_TASK_ID}: Processing ${read}"

# Align reads with minimap2
minimap2 -ax map-ont ${REF} "$read" > "${base}.sam"
echo "Alignment for ${read} done"

# Convert SAM to BAM and sort it
samtools view -@ 4 -bS "${base}.sam" | samtools sort -@ 4 -o "${base}.sorted.bam"
echo "BAM file for ${read} ready"

# Index the sorted BAM file
samtools index "${base}.sorted.bam"
echo "Index for ${read} done"

# Remove SAM file to save space
rm "${base}.sam"

echo "Task ${SLURM_ARRAY_TASK_ID}: Finished ${base}"
```

---

## Step 2: Run Clair3

### Input Data

Mv the minimap output files to the directory of your choice to separate them from the reads. That looks more organised.

Prepare a txt file called `list_bam_all.txt` which lists all indexed BAM files:

| Filename |
|----------|
| SGV002.sorted.bam |
| SGV003.sorted.bam |
| SGV004.sorted.bam |

I downloaded model `r1041_e82_400bps_sup_v420` from the Clair3 website and uploaded to NeSI.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=clair3_calling
#SBATCH --output=slurmlog/%x.%A_%a.out
#SBATCH --error=slurmlog/%x.%A_%a.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00:00
#SBATCH --array=1-57

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis/sam_bam

module purge

# Pull container only if it doesn't exist (run in terminal first time):
# apptainer pull clair3.sif docker://hkubal/clair3:v1.2.0

REF=/nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis/sam_bam/SBV005_medaka.fasta
MODEL_NAME=/nesi/nobackup/massey03742/Salmonella_SNP_calling/r1041_e82_400bps_sup_v420


#If the reference is not indexed yet, run this part first outside of the script and then run the clair script
#module load SAMtools/1.22-GCC-12.3.0
# Index reference only if needed
#[[ -f ${REF}.fai ]] || samtools faidx ${REF}


# Get the BAM file for this array task
bam_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list_bam_all.txt)
base=$(basename "$bam_file" .sorted.bam)

INPUT_BAM=./${bam_file}
OUTPUT_DIR=./clair_out_${base}

apptainer exec \
    clair3.sif \
    /opt/bin/run_clair3.sh \
        --bam_fn=${INPUT_BAM} \
        --ref_fn=${REF} \
        --threads=$SLURM_CPUS_PER_TASK \
        --platform=ont \
        --model_path=${MODEL_NAME} \
        --output=${OUTPUT_DIR} \
        --include_all_ctgs \
        --haploid_precise \
        --no_phasing_for_fa

echo "SNP calling for ${bam_file} done"
```

### Output files inside each folder

| Filename |
|----------|
| full_alignment.vcf.gz |
| merge_output.vcf.gz |
| pileup.vcf.gz |
| run_clair3.log |

---

## Step 3: Extract, rename and filter VCF files

### 3a. Extract and rename VCF files

Because each VCF file is located inside a parent folder with a generic name, we need to extract and rename them according to the isolate name.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=vcf_extract_rename
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis/Clair3_individ

# Set the destination folder
DEST_FOLDER="./vcf_processing"
mkdir -p "$DEST_FOLDER"

# Loop over all clair output folders
for dir in clair_out_*; do
    if [[ -d "$dir" ]]; then
        # Extract the isolate name (e.g., SGV002)
        isolate=$(echo "$dir" | sed 's/clair_out_//')

        # Define the full path to the VCF file
        vcf_file="$dir/merge_output.vcf.gz"
Filter parameters explained:

| Parameter | Explanation |
|-----------|-------------|

```bash
#!/bin/bash -e
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis/vcf_processing

module purge
module load BCFtools/1.22-GCC-12.3.0
module load SAMtools/1.22-GCC-12.3.0

REF=./SBV005_medaka.fasta

# Index reference if needed
if [[ ! -f "${REF}.fai" ]]; then
    echo "Indexing reference..."
    samtools faidx ${REF}
    echo "Reference indexed"
fi

for vcf in *.vcf.gz; do
    # Skip if already processed
    [[ "$vcf" == *".filtered.vcf.gz" ]] && continue

    echo "Processing $vcf"
    sample=$(basename "$vcf" .vcf.gz)

    # Index if needed
    [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]] && bcftools index "$vcf"

    # Filter: QUAL>20, DP>10, SNPs only
    bcftools filter -i 'QUAL>20 && DP>10 && TYPE="snp"' "$vcf" | \
    bcftools norm -f ${REF} -Oz -o "${sample}.filtered.vcf.gz"

    bcftools index "${sample}.filtered.vcf.gz"

    echo "Filtered: ${sample}.filtered.vcf.gz"
done

echo "Finished filtering all VCFs"

```

---

## Step 4: Generate callable position masks

**This is the critical step that was missing from the original pipeline.**

For each isolate, we create a BED file of positions with sufficient coverage (≥10x). Positions below this threshold will be masked as N in the consensus, rather than being incorrectly called as reference.

This mimics what Snippy does internally — only positions with confident calls are included.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=callable_masks
#SBATCH --output=slurmlog/%x.%A_%a.out
#SBATCH --error=slurmlog/%x.%A_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --array=1-366

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis

module purge
module load SAMtools/1.22-GCC-12.3.0

MIN_DEPTH=10
BAM_DIR=./sam_bam
mkdir -p callable_masks

# Get the BAM file for this array task
bam=$(ls ${BAM_DIR}/*.sorted.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(basename "$bam" .sorted.bam)

echo "Task ${SLURM_ARRAY_TASK_ID}: Generating callable mask for $sample..."

# Generate depth file and create BED of callable positions
samtools depth -a "$bam" | \
awk -v min_dp=${MIN_DEPTH} '$3 >= min_dp {print $1"\t"$2-1"\t"$2}' > "callable_masks/${sample}.callable.bed"

# Report coverage stats
total_pos=$(samtools depth -a "$bam" | wc -l)
callable_pos=$(wc -l < "callable_masks/${sample}.callable.bed")
pct=$(echo "scale=2; $callable_pos * 100 / $total_pos" | bc)

echo "$sample: $callable_pos / $total_pos positions callable (${pct}%)"
echo "Task ${SLURM_ARRAY_TASK_ID}: Done"

```

---

## Step 5: Generate consensus FASTA per isolate

Now we generate a consensus sequence for each isolate by:
1. Applying SNPs from the filtered VCF
2. Masking low-coverage positions as N using the callable BED mask

This produces sequences equivalent to Snippy's `.aligned.fa` output.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=generate_consensus
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis/

module purge
module load BCFtools/1.22-GCC-12.3.0
module load BEDTools/2.31.1-GCC-12.3.0

REF=SBV005_medaka.fasta
VCF_DIR=./vcf_processing
MASK_DIR=./callable_masks

mkdir -p consensus_fasta

# Get reference length for creating the inverse mask
ref_length=$(awk '/^>/{if(l) print l; l=0; next}{l+=length}END{print l}' ${REF})
ref_name=$(grep "^>" ${REF} | head -1 | sed 's/>//' | awk '{print $1}')

echo "Reference: $ref_name, Length: $ref_length"

for vcf in ${VCF_DIR}/*.filtered.vcf.gz; do
    sample=$(basename "$vcf" .filtered.vcf.gz)
    callable_bed="${MASK_DIR}/${sample}.callable.bed"

    if [[ ! -f "$callable_bed" ]]; then
        echo "WARNING: No callable mask for $sample, skipping"
        continue
    fi

    echo "Processing $sample..."

    # Create a BED file of the full reference
    echo -e "${ref_name}\t0\t${ref_length}" > temp_full_ref.bed

    # Invert the callable mask to get non-callable (masked) regions
    bedtools subtract -a temp_full_ref.bed -b "$callable_bed" > "temp_${sample}.mask.bed"

    # Generate consensus:
    # -m applies mask (these positions become N)
    # -M N means use N for masked positions
    # -H 1 uses first allele for heterozygous calls (shouldn't happen in haploid)
    bcftools consensus \
        -f ${REF} \
        -m "temp_${sample}.mask.bed" \
        -M N \
        -H 1 \
        "$vcf" | \
    sed "s/^>.*/>$sample/" > "consensus_fasta/${sample}.consensus.fasta"

    # Clean up temp files
    rm "temp_${sample}.mask.bed"

    echo "Generated: consensus_fasta/${sample}.consensus.fasta"
done

rm -f temp_full_ref.bed

# Add reference to the consensus folder (renamed)
sed 's/^>.*/>Reference/' ${REF} > consensus_fasta/Reference.consensus.fasta

echo "All consensus sequences generated"
```

---

## Step 6: Extract core genome alignment using snp-sites

Now we combine all consensus FASTAs and use snp-sites to extract core positions (positions present in all isolates, excluding Ns).

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=snp_sites_core
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Give/Clair/

module purge
module load snp-sites/2.5.1-GCC-11.3.0

mkdir -p alignment_output
mkdir -p alignment_stats

# Concatenate all consensus FASTAs into one file
cat consensus_fasta/*.consensus.fasta > alignment_output/all_isolates.fasta

echo "=== Checking input alignment ==="
echo "Number of sequences:"
grep -c "^>" alignment_output/all_isolates.fasta | tee alignment_stats/input_summary.txt

# Check sequence lengths
module purge
module load SeqKit/2.12.0
echo "Sequence lengths:" | tee -a alignment_stats/input_summary.txt
seqkit fx2tab -nl alignment_output/all_isolates.fasta | tee alignment_stats/input_lengths.tsv

module purge
module load snp-sites/2.5.1-GCC-11.3.0

echo "=== Generating core.full.aln (all core positions, including invariant) ==="
# -c outputs only columns with no gaps/Ns in any sequence (core positions)
snp-sites -c -o alignment_output/core.full.aln alignment_output/all_isolates.fasta

echo "=== Generating core.aln (SNP-only alignment) ==="
# Default output is SNP-only
snp-sites -o alignment_output/core.aln alignment_output/all_isolates.fasta

echo "=== Generating VCF of core SNPs ==="
snp-sites -v -o alignment_output/core.vcf alignment_output/all_isolates.fasta

echo "=== Alignment statistics ==="
module purge
module load SeqKit/2.12.0

echo "Core full alignment (invariant + variant sites):" | tee alignment_stats/core_summary.txt
seqkit stats alignment_output/core.full.aln | tee -a alignment_stats/core_summary.txt

echo "" | tee -a alignment_stats/core_summary.txt
echo "Core SNP alignment (variant sites only):" | tee -a alignment_stats/core_summary.txt
seqkit stats alignment_output/core.aln | tee -a alignment_stats/core_summary.txt

echo "" | tee -a alignment_stats/core_summary.txt
echo "Number of core SNPs:" | tee -a alignment_stats/core_summary.txt
grep -v "^>" alignment_output/core.aln | head -1 | wc -c | tee -a alignment_stats/core_summary.txt

echo "=== Done ==="
echo "Outputs:"
echo "  - alignment_output/core.full.aln  (for Gubbins/BEAST)"
echo "  - alignment_output/core.aln       (SNP-only, for quick trees)"
echo "  - alignment_output/core.vcf       (VCF format)"

```

---

## Step 7: Filter recombination with Gubbins

Run Gubbins on the core full alignment to identify and mask recombinant regions.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=gubbins
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Bovis/Clair3_Bovis/alignment_output

module purge
module load Gubbins/3.2.2-gimkl-2022a-Python-3.10.5
module load IQ-TREE/2.2.2.2-gimpi-2022a

INPUT_ALN=all_isolates.fasta
PREFIX=gubbins_out

# Run Gubbins
run_gubbins.py -p ${PREFIX} ${INPUT_ALN}

echo "=== Finished Gubbins analysis ==="

# Mask recombination sites in the alignment
# This turns recombinant positions to Ns
gff=${PREFIX}.recombination_predictions.gff

mask_gubbins_aln.py --aln ${INPUT_ALN} --gff ${gff} --out ${PREFIX}.masked.fasta

echo "=== Finished masking recombination ==="
echo "Output: ${PREFIX}.masked.fasta (use this for BEAST/TempEst)"
```

---

## Step 8: Generate SNP distance matrix

Use the Gubbins-masked alignment to generate SNP distances.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=snp_dists
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Give/all/alignment_output

module purge
module load nullarbor/2.0.20191013

# Generate SNP distance matrix from masked alignment
snp-dists gubbins_out.masked.fasta > snp_distance_matrix.tsv

echo "SNP distance matrix saved to: snp_distance_matrix.tsv"

# Also generate from unmasked for comparison
snp-dists core.full.aln > snp_distance_matrix_unmasked.tsv

echo "Unmasked SNP distance matrix saved to: snp_distance_matrix_unmasked.tsv"
```

---

## Step 9: Generate phylogenetic tree for TempEst

Generate a maximum likelihood tree from the Gubbins-masked alignment for TempEst analysis.

```bash
#!/bin/bash -e
#SBATCH --account=massey03742
#SBATCH --job-name=iqtree
#SBATCH --output=slurmlog/%x.%j.out
#SBATCH --error=slurmlog/%x.%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

cd /nesi/nobackup/massey03742/Salmonella_SNP_calling/Give/all/alignment_output

module purge
module load IQ-TREE/2.2.2.2-gimpi-2022a

# Run IQ-TREE with model finder and bootstrap
iqtree2 -s gubbins_out.masked.fasta \
    -m MFP \
    -bb 1000 \
    -nt AUTO \
    -pre iqtree_gubbins_masked

echo "=== IQ-TREE finished ==="
echo "Tree file: iqtree_gubbins_masked.treefile (use this in TempEst)"
```

---

## Output Summary

| File | Description | Use for |
|------|-------------|---------|
| `core.full.aln` | Full core alignment (invariant + variant) | Gubbins, BEAST |
| `core.aln` | SNP-only alignment | Quick trees, distance matrices |
| `gubbins_out.masked.fasta` | Recombination-filtered alignment | TempEst, BEAST |
| `snp_distance_matrix.tsv` | Pairwise SNP distances | Epidemiology |
| `iqtree_gubbins_masked.treefile` | ML phylogenetic tree | TempEst |

---


## Link to thesis

- This analysis is linked to section [insert section] of my thesis.


---

*Last updated: 06 February 2026*
