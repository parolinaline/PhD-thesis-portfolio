#Description

Here I explain how I extracted the amino acid alleles from the *cmeRABC* operon in *Campylobacter coli* and *Campylobacter jejuni* genomes and assigned a number to each allele to prepare a amino acid allele profile of the genomes.

## STEP 01
### Annotating the genomes

The script below uses Prokka to annotate the genomes.

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

# Load Prokka module
module purge
module load prokka/1.14.6-apptainer

# Set directories
INPUT_DIR="/home/acalarga/00_nesi_projects/massey03742/sequences_SACNZ/allsequences/jejuni"
OUTPUT_DIR="/home/acalarga/00_nesi_projects/massey03742/sequences_SACNZ/allsequences/jejuni/prokka_output"

mkdir -p ${OUTPUT_DIR}

# Each task processes 10 genomes
BATCH_SIZE=10
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( START + BATCH_SIZE - 1 ))

# Get the list of fasta files
FILELIST=($(ls ${INPUT_DIR}/*.fasta))
TOTAL=${#FILELIST[@]}

for i in $(seq ${START} ${END}); do
    # Stop if we've exceeded the total number of files
    INDEX=$((i - 1))
    if [ ${INDEX} -ge ${TOTAL} ]; then
        echo "No more files to process, exiting."
        break
    fi

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

echo "Task ${SLURM_ARRAY_TASK_ID} complete."
```

## STEP 02
### BLAST the genomes annotations against the cmeABCR reference

The script below will extract all annotated genes from the Prokka output, rename the headers so each one has the isolate's name which is the parent folder name.
Then, I made the BLASTp database using the output file with all proteins from all annotations.
Get the genes cmeA, cmeB, cmeC, cmeR from the *Campylobacter jejuni* NCTC11168 for the *C. jejuni* genomes.
For the *Campylobacter coli* genomes I used the reference genome *Campylobacter coli* isolate SC0001 from the hybrid assembly I obtained.

#### *Campylobacter jejuni* cmeABCR reference.

I saved as cmeABCR.faa in *C. jejuni*`s folder:

```bash
>CmeA
MKLFQKNTILALGVVLLLTACSKEEAPKIQMPPQPVTTMSAKSEDLPLSFTYPAKLVSDYDVIIKPQVSGVIENKLFKAGDKVKKGQTLFIIEQDKFKASVDSAYGQALMAKATFENASKDFNRSKALFSKSAISQKEYDSSLATFNNSKASLASARAQLANARIDLDHTEIKAPFDGTIGDALVNIGDYVSASTTELVRVTNLNPIYADFFISDTDKLNLVRNTQNGKWDLDSIHANLNLNGETVQGKLYFIDSVIDANSGTVKAKAIFDNNNSTLLPGAFATITSEGFIQKNGFKVPQIAVKQNQNDVYVLLVKNGKVEKSSVHISYQNNEYAIIDKGLQNGDKIILDNFKKIQVGSEVKEIGAQ
>CmeB
MFSKFFIERPVFASVVAIIISLAGAIGLTNLPIEQYPSLTPPTVKVSATYTGADAQTIASTVASPIEDAINGADNMIYMDSTSSSSGTMSLTVYFDIGTDPDQATIDVNNRISAATAKMPDAVKKLGVTVRKTSSTTLAAISMYSSDGSMSAVDVYNYITLNVLDELKRVPGVGDANAIGNRNYSLRIWLKPDLLNKFGITATDVISAVNDQNAQYATGKIGEEPVTQKSPYVYSITMQGRLQNPSEFENIILRTNNDGSFLRLKDVADVEIGSQQYSSQGRLNGNDAVPIMINLQSGANALHTAELVQAKMQELSKNFPKGLTYKIPYDTTKFVIESIKEVVKTFVEALILVIIVMYMFLKNFRATLIPMIAVPVSLLGTFAGLYVLGFSINLLTLFALILAIGIVVDDAIIVVENIDRILHENEQISVKDAAIQAMQEVSSPVISIVLVLCAVFVPVSFISGFVGEIQRQFALTLAISVTISGFVALTLTPSLCALFLRRNEGEPFKFVKKFNDFFDWSTSVFSAGVAYILKRTIRFVLIFCIMLGAIFYLYKAVPSSLVPEEDQGLMIGIINLPSASALHRTISEVDHISQEVLKTNGVKDAMAMIGFDLFTSSLKENAAAMFIGLKDWKDRNVSADEIAMELNKKFAFDRNASSIFIGLPPIPGLSITGGFEMYVQNKSGKSYDEIQKDVNKLVAVANQRKELSRVRTTLDTTFPQYKLIIDRDKLKHYNLNMQDVFNTMNATIGTYYVNDFSMLGKNFQVNIRAKGDFRNTQDALKNIFVRSNDGKMIPLDSFLTLQRSSGPDDVKRFNLFPAAQVQGQPAPGYTSGQAIEAIAQVAKETLGDDYSIAWSGSAYQEVSSKGTASYAFALGMIFVFLILAAQYERWLIPLAVVTAVPFAVFGSFLLVYLRGFSNDIYFQTGLLLLIGLSAKNAILIVEFAMEERFKKGKGVFEAAVAAAKLRFRPIIMTSLAFTFGVLPMIFATGAGSASRHSLGTGLIGGMIAASTLAIFFVPLFFYLLENFNEWLDKKRGKVHE
>CmeC
MNKIISISAIASFTLLISACSLSPNLNIPEANYSIDNKLGALSWEKENNSSITKNWWKDFDDENLNKVVDLALKNNNDLKLAFIHMEQAAAQLGIDFSSLLPKFDGSASGSRAKTAINAPSNRTGEVSYGNDFKMGLNLSYEIDLWGKYRDTYRASKSGFKASEYDYEAARLSVISNTVQTYFNLVNAYENENALKEAYKSAKEIYRINDEKFQVGAVGEYELAQARANLESMALQYNEAKLNKENYLKALKILTSNDLNDILYKNQSYQVFNLKEFDIPTGISSTILLQRPDIGSSLEKLTQQNYLVGVARTAFLPSLSLTGLLGFESGDLDTLVKGGSKTWNIGGNFTLPIFHWGEIYQNVNLAKLNKDEAFVNYQNTLITAFGEIRYALVARKTIRLQYDNAQASEQSYKRIYEIAKERYDIGEMSLQDYLEARQNWLNAAVAFNNIKYSYANSIVDVIKAFGGGFEQSEDTSKNIKEESKNLDMSFRE
>CmeR
MNSNRTPSQKVLARQEKIKAVALELFLTKGYQETSLSDIIKLSG GSYSNIYDGFKSKEGLFFEILDDICKKHFHLIYSKTQEIENGTLKEILTSFGLAFIE  FNQPEAVAFGKIIYSQVYDKDRHLANWIENNQQNFSYNILMGFFKQQNNSYMKKNA   LAVLFCTMLKEPYHHLNVLINAPLKNKKEQKEHVEFVVNVFLNGINSSKA
```
#### *Campylobacter coli* cmeABCR reference.


I saved as cmeABCR_coli.faa in *C. coli*`s folder:
```bash
>CmeA
MNLFQKNTLLLLSALFLFSACSKEEAPQKQTPPQSVSTMSAKAENLPLNFTYPAKLVSDYDVIIKPQVSGVIVEKLFKAGDLIKKGQTLFIIEQDKFKASVNSAYGKALMARANFDNASKDYNRSKTLYNKGAISQKEYDSALANFNNTKANLTSARADLENARIDLAYTEIKAPFDGIVGDALINIGDYVSSSSTELVRITNLNPIYADFYISDTDKLNIVRNTQDGKWDLSNIYADLNLNGEVVKGKLYFIDSVIDANSGTVKAKAIFDNNDSTLLPGAFATITSNGFIQKNGFKIPQIAIKQDQNEVYVFLLKEGKVAKAPVHISYQDNEYAIIDKGLQNGDKIILDNFKKIRLGSEVKEVGAQ
>CmeB
MFSKFFIERPVFASVVAIIISLAGAIGLVNLPIEQYPSLTPPTVKVSATYTGADAQTIATTVANPIEDAINGADNMIYMDSTSSSSGTMNLTVYFNIGTDPDQATIDVNNRISAATAKMPDEVKKLGVTVRKTSSTTLAAISMYSSNGSMSAVDVYNYITLNILDELKRVPGVGDANAMGNRNYSLRIWLKPDLLNKFGATATDVIAAVNDQNAQYATGKIGEEPVTQKSPYVYSITMQGRLQSPSEFENIILRTNEDGSFLRLKDVADVEIGSQQYKSQGRLNGNDAVPIMINLQSGANALNTAKLVEAKMQELSKSFPEGLEYKIPYDTTKFVIESIKEVIKTFVEALILVIIVMYMFLKNFRSTLIPMIAVPVSLLGTFAGLYLLGFSINLLTLFALILAIGIVVDDAIIVVENIDRILHEDEKISVKDAAIQAMQEVSSPVISIVLVLCAVFIPVSFISGFVGEIQKQFALTLAISVTISGFVALTLTPSLCALFLRRNESKPFYIVQKFNDFFDWSTSIFSAGVAYMLKRTIRFVLIFCIMLGAIFYLYKQVPGSLVPEEDQGLMIGIVNLPSASALHRTISEVDSMSQEILKTNGVKDAMAMIGFDLFTSSLKENAAAMFIGLEDWKDRNVSADEIIMELNKKFAPDRNAASVFRGLPPIPGLSITGGFEMYVQNKSGKSYDQIQEDVNKLVAAANQREELYGVRTTLDTSFPQYKLIIDRDKLKHFNLNMQDVFSTMNATIGTYYVNDFTMLGKNFQVNIRAKGDFINTQNALKNIFVRSNDGKMIPLDSFLTLQRSSGPDDVKRFNIFPAAQIQGQPAPGYTSGQAIEAISQVAQETLSDDYSIAWSGSAYQEVSSKGTGSYAFALGMVFVFLILAAQYERWLIPLAVVTAVPFAVFGSFLLVYLRGFTNDIYFQTGLLLLIGLSAKNAILIVEFAMEERFKKGKKIFDAAIEAAKLRFRPIVMTSLAFTFGVLPMIFATGAGSASRHSLGTGLIGGMIAASTLAIFFVPLFFYLLENFNEWLDKKRGKVHE
>CmeC
MNKIISISVIASLSLFISACSLSPNLNIPEANYSLDNKLGALSWEKENNHTLSKEWWKEFDDEKLNHVVDLALQNNNDLKLAFIHMEQAAAQLGIDFSDLLPKLDGSASANRAKTPTNAPSNRTGGVSYGNDFKMGLDLSYEIDLWGRYRDTYRASRSGFKASQYDYEAARLSIVSSVVQTYFNFVNANENEKALKDAYESAQEIYQINYEKFQVGAVGEYEIAQSRANLESTALQYNEAKLNKENYLKALKILTSNDLNDILYQDQAYQVFKLKDFDIPSGISSTILLQRPDIGSSLEKLTQQNYLVGVARTAFLPNLSLTGLLGFESGDLNTLVEAGSRTWSIGGNFVMPIFHWGEIYQNVNLAKLNKDEAFVNYQNTLVAAFGEIRYALIARKTIRLQYNNAQASELSYKRIYEISKERYDVGEMSLQDYLQARQDWLNATVAFNNTKYSYANSIVNVIKAFGGGFEQGKNISKNIEEESKTLDMSFRE
>CmeR
MNPNKTPSKKVLARREKIKNVAFDLFLTKGFQETSLSDIIKLSGGSYSNIYDSFNSKEGLFFEILDDVCKKHFDLIASQTQTIKDKNLKEFLTSFGLTFVDIFNQVQTVAFGKIIFSQVYDKHKHVENWIENNQKIFSYNILIELFKKQDSSYISNNAQKLAMLFCAMLREPYHSLNVLADTPLMNKQEQKEHVEFIVNIFLKGIENKTSI

```

### Run the following script:

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

# Concatenate all Prokka .faa files with genome name in header
for dir in /nesi/project/massey03742/sequences_SACNZ/allsequences/jejuni/prokka_output/SC*/; do
    genome=$(basename "$dir")
    sed "s/^>/>$genome|/" "$dir"*.faa
done > all_proteins.faa

# Make BLAST db
makeblastdb -in all_proteins.faa -dbtype prot -out all_proteins_db

# BLAST your references
blastp -query cmeABCR.faa -db all_proteins_db \
  -outfmt "6 qseqid sseqid pident qcovs length sseq" \
  -evalue 1e-20 -max_target_seqs 5000 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -out cme_hits.tsv
```

### Filtering

Because cme genes have ortholog genes within *Campylobacter* genomes, cmeB and cmeC have two hits in the genome, one with low % of identity and the other one with >80% which is the real hit. So, we are going to filter the output to only save the hits with percentage of identity >80%:

```bash
awk -F'\t' '$3 >= 80' cme_hits.tsv > cme_hits_filtered.tsv
```
## STEP03
### Assigning amino acid alleles to the genomes

The following script will get the BLASTp results and assign an ID to each allele.

```python
#!/usr/bin/env python3
"""
Extract CmeABCR from BLAST results against Prokka proteins,
assign AA allele IDs, and generate a profile table.
"""

import csv
from collections import defaultdict

blast_file = "cme_hits_filtered.tsv"
output_table = "cmeABCR_allele_profiles.tsv"
output_fastas = "cme_allele_fastas"

# Thresholds
MIN_PIDENT = 70   # minimum % identity
MIN_QCOV = 80     # minimum query coverage %

target_genes = ["CmeA", "CmeB", "CmeC", "CmeR"]

import os
os.makedirs(output_fastas, exist_ok=True)

# gene -> {aa_seq: allele_id}
allele_dicts = {g: {} for g in target_genes}
allele_counters = {g: 1 for g in target_genes}
# genome -> {gene: allele_id}
profiles = defaultdict(dict)
# track best hit per genome per gene (by pident)
best_hits = defaultdict(dict)  # (genome, gene) -> (pident, aa_seq)

with open(blast_file) as f:
    for line in f:
        cols = line.strip().split("\t")
        gene = cols[0]       # e.g., CmeA
        subject = cols[1]    # e.g., genome_name|locus_tag
        pident = float(cols[2])
        qcov = float(cols[3])
        aa_seq = cols[5].replace("-", "").rstrip("*")

        if pident < MIN_PIDENT or qcov < MIN_QCOV:
            continue

        genome = subject.split("|")[0]
        key = (genome, gene)

        # Keep only best hit per genome per gene
        if key not in best_hits or pident > best_hits[key][0]:
            best_hits[key] = (pident, aa_seq)

# Assign allele IDs
for (genome, gene), (pident, aa_seq) in sorted(best_hits.items()):
    if aa_seq not in allele_dicts[gene]:
        allele_dicts[gene][aa_seq] = allele_counters[gene]
        allele_counters[gene] += 1
    profiles[genome][gene] = allele_dicts[gene][aa_seq]

# Write profile table
with open(output_table, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["genome"] + target_genes)
    for genome in sorted(profiles):
        row = [genome] + [profiles[genome].get(g, "-") for g in target_genes]
        writer.writerow(row)

# Write allele FASTAs
for gene in target_genes:
    with open(os.path.join(output_fastas, f"{gene}_alleles.faa"), "w") as f:
        for aa_seq, aid in sorted(allele_dicts[gene].items(), key=lambda x: x[1]):
            f.write(f">{gene}_allele_{aid}\n{aa_seq}\n")
    print(f"{gene}: {len(allele_dicts[gene])} unique alleles")

print(f"\nProfiles for {len(profiles)} genomes written to {output_table}")
```
