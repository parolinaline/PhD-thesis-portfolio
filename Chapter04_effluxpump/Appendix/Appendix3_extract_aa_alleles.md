# Amino Acid Allele Extraction from the *cmeRABC* Operon

This document describes the pipeline used to extract amino acid sequences from the *cmeR*, *cmeA*, *cmeB*, and *cmeC* genes in *Campylobacter jejuni* and *Campylobacter coli* genomes, assign a unique numeric ID to each distinct allele, and generate per-genome allele profiles.

---

## Overview

The pipeline consists of three main steps:

1. **Genome annotation** — annotate all genomes using Prokka to obtain predicted protein sequences (`.faa` files).
2. **BLASTp search** — query species-specific CmeABCR reference sequences against a combined database of all annotated proteins.
3. **Allele assignment** — parse BLAST hits, deduplicate by best hit per genome per gene, assign a numeric ID to each unique amino acid sequence, and write the allele profiles.

```
Genome assemblies (.fasta)
        │
        ▼
   [STEP 01] Prokka annotation
        │
        ▼ per-genome .faa files
   [STEP 02] Concatenate proteins → BLASTp vs. CmeABCR references
        │
        ▼ cme_hits_filtered.tsv (pident ≥ 80%)
   [STEP 03] Python: allele assignment → cmeABCR_allele_profiles.tsv
                                       → cme_allele_fastas/
```

---

## Reference sequences

Species-specific reference sequences for CmeA, CmeB, CmeC, and CmeR were used as BLAST queries.

- **_C. jejuni_**: reference genes retrieved from *C. jejuni* NCTC11168, saved as `cmeABCR.faa`.
- **_C. coli_**: reference genes retrieved from the hybrid-assembly isolate SC0001, saved as `cmeABCR_coli.faa`.

<details>
<summary><em>C. jejuni</em> reference sequences — <code>cmeABCR.faa</code></summary>

```fasta
>CmeA
MKLFQKNTILALGVVLLLTACSKEEAPKIQMPPQPVTTMSAKSEDLPLSFTYPAKLVSDYDVIIKPQVSGVIENKLFKAGDKVKKGQTLFIIEQDKFKASVDSAYGQALMAKATFENASKDFNRSKALFSKSAISQKEYDSSLATFNNSKASLASARAQLANARIDLDHTEIKAPFDGTIGDALVNIGDYVSASTTELVRVTNLNPIYADFFISDTDKLNLVRNTQNGKWDLDSIHANLNLNGETVQGKLYFIDSVIDANSGTVKAKAIFDNNNSTLLPGAFATITSEGFIQKNGFKVPQIAVKQNQNDVYVLLVKNGKVEKSSVHISYQNNEYAIIDKGLQNGDKIILDNFKKIQVGSEVKEIGAQ
>CmeB
MFSKFFIERPVFASVVAIIISLAGAIGLTNLPIEQYPSLTPPTVKVSATYTGADAQTIASTVASPIEDAINGADNMIYMDSTSSSSGTMSLTVYFDIGTDPDQATIDVNNRISAATAKMPDAVKKLGVTVRKTSSTTLAAISMYSSDGSMSAVDVYNYITLNVLDELKRVPGVGDANAIGNRNYSLRIWLKPDLLNKFGITATDVISAVNDQNAQYATGKIGEEPVTQKSPYVYSITMQGRLQNPSEFENIILRTNNDGSFLRLKDVADVEIGSQQYSSQGRLNGNDAVPIMINLQSGANALHTAELVQAKMQELSKNFPKGLTYKIPYDTTKFVIESIKEVVKTFVEALILVIIVMYMFLKNFRATLIPMIAVPVSLLGTFAGLYVLGFSINLLTLFALILAIGIVVDDAIIVVENIDRILHENEQISVKDAAIQAMQEVSSPVISIVLVLCAVFVPVSFISGFVGEIQRQFALTLAISVTISGFVALTLTPSLCALFLRRNEGEPFKFVKKFNDFFDWSTSVFSAGVAYILKRTIRFVLIFCIMLGAIFYLYKAVPSSLVPEEDQGLMIGIINLPSASALHRTISEVDHISQEVLKTNGVKDAMAMIGFDLFTSSLKENAAAMFIGLKDWKDRNVSADEIAMELNKKFAFDRNASSIFIGLPPIPGLSITGGFEMYVQNKSGKSYDEIQKDVNKLVAVANQRKELSRVRTTLDTTFPQYKLIIDRDKLKHYNLNMQDVFNTMNATIGTYYVNDFSMLGKNFQVNIRAKGDFRNTQDALKNIFVRSNDGKMIPLDSFLTLQRSSGPDDVKRFNLFPAAQVQGQPAPGYTSGQAIEAIAQVAKETLGDDYSIAWSGSAYQEVSSKGTASYAFALGMIFVFLILAAQYERWLIPLAVVTAVPFAVFGSFLLVYLRGFSNDIYFQTGLLLLIGLSAKNAILIVEFAMEERFKKGKGVFEAAVAAAKLRFRPIIMTSLAFTFGVLPMIFATGAGSASRHSLGTGLIGGMIAASTLAIFFVPLFFYLLENFNEWLDKKRGKVHE
>CmeC
MNKIISISAIASFTLLISACSLSPNLNIPEANYSIDNKLGALSWEKENNSSITKNWWKDFDDENLNKVVDLALKNNNDLKLAFIHMEQAAAQLGIDFSSLLPKFDGSASGSRAKTAINAPSNRTGEVSYGNDFKMGLNLSYEIDLWGKYRDTYRASKSGFKASEYDYEAARLSVISNTVQTYFNLVNAYENENALKEAYKSAKEIYRINDEKFQVGAVGEYELAQARANLESMALQYNEAKLNKENYLKALKILTSNDLNDILYKNQSYQVFNLKEFDIPTGISSTILLQRPDIGSSLEKLTQQNYLVGVARTAFLPSLSLTGLLGFESGDLDTLVKGGSKTWNIGGNFTLPIFHWGEIYQNVNLAKLNKDEAFVNYQNTLITAFGEIRYALVARKTIRLQYDNAQASEQSYKRIYEIAKERYDIGEMSLQDYLEARQNWLNAAVAFNNIKYSYANSIVDVIKAFGGGFEQSEDTSKNIKEESKNLDMSFRE
>CmeR
MNSNRTPSQKVLARQEKIKAVALELFLTKGYQETSLSDIIKLSGGSYSNIYDGFKSKEGLFFEILDDICKKHFHLIYSKTQEIENGTLKEILTSFGLAFIE
FNQPEAVAFGKIIYSQVYDKDRHLANWIENNQQNFSYNILMGFFKQQNNSYMKKNA
LAVLFCTMLKEPYHHLNVLINAPLKNKKEQKEHVEFVVNVFLNGINSSKA
```
</details>

<details>
<summary><em>C. coli</em> reference sequences — <code>cmeABCR_coli.faa</code></summary>

```fasta
>CmeA
MNLFQKNTLLLLSALFLFSACSKEEAPQKQTPPQSVSTMSAKAENLPLNFTYPAKLVSDYDVIIKPQVSGVIVEKLFKAGDLIKKGQTLFIIEQDKFKASVNSAYGKALMARANFDNASKDYNRSKTLYNKGAISQKEYDSALANFNNTKANLTSARADLENARIDLAYTEIKAPFDGIVGDALINIGDYVSSSSTELVRITNLNPIYADFYISDTDKLNIVRNTQDGKWDLSNIYADLNLNGEVVKGKLYFIDSVIDANSGTVKAKAIFDNNDSTLLPGAFATITSNGFIQKNGFKIPQIAIKQDQNEVYVFLLKEGKVAKAPVHISYQDNEYAIIDKGLQNGDKIILDNFKKIRLGSEVKEVGAQ
>CmeB
MFSKFFIERPVFASVVAIIISLAGAIGLVNLPIEQYPSLTPPTVKVSATYTGADAQTIATTVANPIEDAINGADNMIYMDSTSSSSGTMNLTVYFNIGTDPDQATIDVNNRISAATAKMPDEVKKLGVTVRKTSSTTLAAISMYSSNGSMSAVDVYNYITLNILDELKRVPGVGDANAMGNRNYSLRIWLKPDLLNKFGATATDVIAAVNDQNAQYATGKIGEEPVTQKSPYVYSITMQGRLQSPSEFENIILRTNEDGSFLRLKDVADVEIGSQQYKSQGRLNGNDAVPIMINLQSGANALNTAKLVEAKMQELSKSFPEGLEYKIPYDTTKFVIESIKEVIKTFVEALILVIIVMYMFLKNFRSTLIPMIAVPVSLLGTFAGLYLLGFSINLLTLFALILAIGIVVDDAIIVVENIDRILHEDEKISVKDAAIQAMQEVSSPVISIVLVLCAVFIPVSFISGFVGEIQKQFALTLAISVTISGFVALTLTPSLCALFLRRNESKPFYIVQKFNDFFDWSTSIFSAGVAYMLKRTIRFVLIFCIMLGAIFYLYKQVPGSLVPEEDQGLMIGIVNLPSASALHRTISEVDSMSQEILKTNGVKDAMAMIGFDLFTSSLKENAAAMFIGLEDWKDRNVSADEIIMELNKKFAPDRNAASVFRGLPPIPGLSITGGFEMYVQNKSGKSYDQIQEDVNKLVAAANQREELYGVRTTLDTSFPQYKLIIDRDKLKHFNLNMQDVFSTMNATIGTYYVNDFTMLGKNFQVNIRAKGDFINTQNALKNIFVRSNDGKMIPLDSFLTLQRSSGPDDVKRFNIFPAAQIQGQPAPGYTSGQAIEAISQVAQETLSDDYSIAWSGSAYQEVSSKGTGSYAFALGMVFVFLILAAQYERWLIPLAVVTAVPFAVFGSFLLVYLRGFTNDIYFQTGLLLLIGLSAKNAILIVEFAMEERFKKGKKIFDAAIEAAKLRFRPIVMTSLAFTFGVLPMIFATGAGSASRHSLGTGLIGGMIAASTLAIFFVPLFFYLLENFNEWLDKKRGKVHE
>CmeC
MNKIISISVIASLSLFISACSLSPNLNIPEANYSLDNKLGALSWEKENNHTLSKEWWKEFDDEKLNHVVDLALQNNNDLKLAFIHMEQAAAQLGIDFSDLLPKLDGSASANRAKTPTNAPSNRTGGVSYGNDFKMGLDLSYEIDLWGRYRDTYRASRSGFKASQYDYEAARLSIVSSVVQTYFNFVNANENEKALKDAYESAQEIYQINYEKFQVGAVGEYEIAQSRANLESTALQYNEAKLNKENYLKALKILTSNDLNDILYQDQAYQVFKLKDFDIPSGISSTILLQRPDIGSSLEKLTQQNYLVGVARTAFLPNLSLTGLLGFESGDLNTLVEAGSRTWSIGGNFVMPIFHWGEIYQNVNLAKLNKDEAFVNYQNTLVAAFGEIRYALIARKTIRLQYNNAQASELSYKRIYEISKERYDVGEMSLQDYLQARQDWLNATVAFNNTKYSYANSIVNVIKAFGGGFEQGKNISKNIEEESKTLDMSFRE
>CmeR
MNPNKTPSKKVLARREKIKNVAFDLFLTKGFQETSLSDIIKLSGGSYSNIYDSFNSKEGLFFEILDDVCKKHFDLIASQTQTIKDKNLKEFLTSFGLTFVDIFNQVQTVAFGKIIFSQVYDKHKHVENWIENNQKIFSYNILIELFKKQDSSYISNNAQKLAMLFCAMLREPYHSLNVLADTPLMNKQEQKEHVEFIVNIFLKGIENKTSI
```
</details>

---

## Step 01 — Genome annotation with Prokka

All genome assemblies (`.fasta`) were annotated using [Prokka v1.14.6](https://github.com/tseemann/prokka) to obtain predicted protein sequences. The job was submitted as a SLURM array on NeSI, processing genomes in batches of 10.

**Script:** `01_prokka_annotation.sh`

```bash
#!/bin/bash -e
#SBATCH --job-name=prokka_cjejuni
#SBATCH --account=massey03742
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-107%50
#SBATCH --output=slurmlog/prokka_%A_%a.out
#SBATCH --error=slurmlog/prokka_%A_%a.err

module purge
module load prokka/1.14.6-apptainer

INPUT_DIR="/home/acalarga/00_nesi_projects/massey03742/sequences_SACNZ/allsequences/jejuni"
OUTPUT_DIR="${INPUT_DIR}/prokka_output"
mkdir -p "${OUTPUT_DIR}"

# Process genomes in batches of 10 per array task
BATCH_SIZE=10
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( START + BATCH_SIZE - 1 ))

FILELIST=($(ls ${INPUT_DIR}/*.fasta))
TOTAL=${#FILELIST[@]}

for i in $(seq ${START} ${END}); do
    INDEX=$((i - 1))
    [ ${INDEX} -ge ${TOTAL} ] && { echo "No more files to process."; break; }

    FASTA=${FILELIST[${INDEX}]}
    SAMPLE=$(basename ${FASTA} .fasta)

    echo "Task ${SLURM_ARRAY_TASK_ID} | [${i}/${TOTAL}] Annotating: ${SAMPLE}"

    prokka \
        --outdir ${OUTPUT_DIR}/${SAMPLE} \
        --prefix ${SAMPLE} \
        --genus Campylobacter \
        --species jejuni \
        --kingdom Bacteria \
        --cpus ${SLURM_CPUS_PER_TASK} \
        --force \
        --addgenes \
        --mincontiglen 200 \
        --rfam \
        ${FASTA}

    echo "Task ${SLURM_ARRAY_TASK_ID} | [${i}/${TOTAL}] Done: ${SAMPLE}"
done
```

**Key Prokka parameters:**

| Parameter | Value | Purpose |
|---|---|---|
| `--genus` / `--species` | *Campylobacter jejuni* | Species-specific annotation database |
| `--addgenes` | — | Include gene names in annotation output |
| `--mincontiglen` | 200 bp | Exclude very short contigs |
| `--rfam` | — | Enable ncRNA annotation via Rfam |
| `--array` | 1–107 (50 at a time) | Parallelise across all genomes |

> **Note:** The same script was used for *C. coli* genomes with `--species coli` and the appropriate input directory.

---

## Step 02 — BLASTp search against CmeABCR references

### Rationale for filtering

CmeB and CmeC each have a paralogue within the *Campylobacter* genome (CmeE and CmeF, respectively), which produces a secondary BLAST hit with low sequence identity. Filtering to hits with **≥ 80% identity** reliably retains only the true *cmeABC* hits.

### Script: `02_blastp_cmeABCR.sh`

```bash
#!/bin/bash -e
#SBATCH --job-name=BLASTp
#SBATCH --account=massey03742
#SBATCH --time=12:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurmlog/%j.out
#SBATCH --error=slurmlog/%j.err

module purge
module load BLAST/2.16.0-GCC-12.3.0

# 1. Concatenate all Prokka protein files, prepending the genome name to each header
for dir in /nesi/project/massey03742/sequences_SACNZ/allsequences/jejuni/prokka_output/SC*/; do
    genome=$(basename "$dir")
    sed "s/^>/>$genome|/" "$dir"*.faa
done > all_proteins.faa

# 2. Build the BLAST protein database
makeblastdb -in all_proteins.faa -dbtype prot -out all_proteins_db

# 3. Query CmeABCR references against the combined database
blastp \
    -query cmeABCR.faa \
    -db all_proteins_db \
    -outfmt "6 qseqid sseqid pident qcovs length sseq" \
    -evalue 1e-20 \
    -max_target_seqs 5000 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out cme_hits.tsv

# 4. Retain only hits with pident >= 80%
awk -F'\t' '$3 >= 80' cme_hits.tsv > cme_hits_filtered.tsv
```

**Output columns in `cme_hits_filtered.tsv`:**

| Column | Field | Description |
|---|---|---|
| 1 | `qseqid` | Query gene name (e.g., `CmeA`) |
| 2 | `sseqid` | Subject ID in format `genome\|locus_tag` |
| 3 | `pident` | Percentage amino acid identity |
| 4 | `qcovs` | Query coverage (%) |
| 5 | `length` | Alignment length |
| 6 | `sseq` | Subject amino acid sequence (aligned region) |

---

## Step 03 — Allele ID assignment

### Logic

For each genome and each gene, the script:

1. Retains the **best BLAST hit** (highest `pident`), applying minimum thresholds of `pident ≥ 70%` and `qcovs ≥ 80%`.
2. Strips gap characters (`-`) and stop codons (`*`) from the subject sequence.
3. Assigns a **numeric allele ID** to each unique amino acid sequence, incrementing the counter per gene.
4. Outputs an allele profile table (one row per genome) and a FASTA file per gene containing all distinct alleles.

> **Note:** The minimum identity threshold in the Python script (`MIN_PIDENT = 70`) is intentionally lower than the Bash pre-filter (`pident ≥ 80`), which has already removed paralogue hits. This provides a safety margin if the pipeline is rerun with a less strict pre-filter.

### Script: `03_assign_alleles.py`

```python
#!/usr/bin/env python3
"""
Parse BLASTp results (cme_hits_filtered.tsv) to assign numeric amino acid
allele IDs to each CmeABCR gene per genome, and write:
  - cmeABCR_allele_profiles.tsv  : per-genome allele profile table
  - cme_allele_fastas/           : per-gene FASTA files of all distinct alleles
"""

import csv
import os
from collections import defaultdict

# --- Parameters ---
blast_file    = "cme_hits_filtered.tsv"
output_table  = "cmeABCR_allele_profiles.tsv"
output_fastas = "cme_allele_fastas"
MIN_PIDENT    = 70   # minimum % amino acid identity
MIN_QCOV      = 80   # minimum query coverage (%)
target_genes  = ["CmeA", "CmeB", "CmeC", "CmeR"]

os.makedirs(output_fastas, exist_ok=True)

# --- Data structures ---
allele_dicts    = {g: {} for g in target_genes}   # gene -> {aa_seq: allele_id}
allele_counters = {g: 1 for g in target_genes}    # gene -> next allele ID
profiles        = defaultdict(dict)               # genome -> {gene: allele_id}
best_hits       = {}                              # (genome, gene) -> (pident, aa_seq)

# --- Pass 1: collect best hit per genome per gene ---
with open(blast_file) as f:
    for line in f:
        cols   = line.strip().split("\t")
        gene   = cols[0]
        subject= cols[1]                          # genome_name|locus_tag
        pident = float(cols[2])
        qcov   = float(cols[3])
        aa_seq = cols[5].replace("-", "").rstrip("*")

        if pident < MIN_PIDENT or qcov < MIN_QCOV:
            continue

        genome = subject.split("|")[0]
        key    = (genome, gene)

        if key not in best_hits or pident > best_hits[key][0]:
            best_hits[key] = (pident, aa_seq)

# --- Pass 2: assign allele IDs ---
for (genome, gene), (pident, aa_seq) in sorted(best_hits.items()):
    if aa_seq not in allele_dicts[gene]:
        allele_dicts[gene][aa_seq] = allele_counters[gene]
        allele_counters[gene] += 1
    profiles[genome][gene] = allele_dicts[gene][aa_seq]

# --- Write profile table ---
with open(output_table, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["genome"] + target_genes)
    for genome in sorted(profiles):
        row = [genome] + [profiles[genome].get(g, "-") for g in target_genes]
        writer.writerow(row)

# --- Write per-gene allele FASTAs ---
for gene in target_genes:
    fasta_path = os.path.join(output_fastas, f"{gene}_alleles.faa")
    with open(fasta_path, "w") as f:
        for aa_seq, aid in sorted(allele_dicts[gene].items(), key=lambda x: x[1]):
            f.write(f">{gene}_allele_{aid}\n{aa_seq}\n")
    print(f"{gene}: {len(allele_dicts[gene])} unique alleles")

print(f"\nProfiles for {len(profiles)} genomes written to {output_table}")
```

### Output files

| File | Description |
|---|---|
| `cmeABCR_allele_profiles.tsv` | Tab-separated table: one row per genome, columns for CmeA / CmeB / CmeC / CmeR allele IDs. Missing genes are coded as `-`. |
| `cme_allele_fastas/CmeA_alleles.faa` | FASTA of all unique CmeA amino acid alleles, named `CmeA_allele_1`, `CmeA_allele_2`, … |
| `cme_allele_fastas/CmeB_alleles.faa` | Same for CmeB |
| `cme_allele_fastas/CmeC_alleles.faa` | Same for CmeC |
| `cme_allele_fastas/CmeR_alleles.faa` | Same for CmeR |

Example of `cmeABCR_allele_profiles.tsv`:

```
genome      CmeA    CmeB    CmeC    CmeR
SC0001      1       1       1       1
SC0002      2       1       1       2
SC0003      1       3       2       1
...
```

---

## Filtering thresholds summary

| Filter | Threshold | Applied at |
|---|---|---|
| E-value | ≤ 1×10⁻²⁰ | BLASTp (`-evalue`) |
| % identity (pre-filter) | ≥ 80% | `awk` post-BLAST |
| % identity (Python) | ≥ 70% | `MIN_PIDENT` |
| Query coverage (Python) | ≥ 80% | `MIN_QCOV` |
| Hit selection | Best `pident` per genome per gene | `best_hits` dict |
