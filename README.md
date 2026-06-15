# LGRASP_Assemblers
Bambu / FLAIR / IsoQuant / StringTie /  
Veuillez trouver les informations associées à l'évaluation des outils d'assemblage présentés dans le mémoire de M2 **"De la détection des isoformes à la classification : Apport du séquençage RNA ciblé long-read au diagnostic du syndrome HBOC"** (section 2.2.2 - page 9) — Appoline NABI

---

# 🇫🇷 SECTION FRANÇAISE

---

## ALIGNEMENT DES FASTQ SUR LE GÉNOME DE RÉFÉRENCE (hg38) — MINIMAP2 en mode splice

```bash
minimap2 \
  -ax splice \
  --MD \
  -t {threads} \
  {params.ref} \
  -o {output.sam} \
  {input}

samtools view \
  -b \
  -@ {threads} \
  -o {output.sorted_bam} \
  {output.sam}

samtools sort \
  -@ {threads} \
  -o {output.bam} \
  {output.sorted_bam}

samtools index \
  -@ {threads} \
  {output.bam}
```

> **Note pour FLAIR** : suivi des recommandations issues de la documentation, suppression des alignements secondaires.

```bash
minimap2 \
  -ax splice \
  -s 80 \
  -G 200k \
  --MD \
  --secondary=no \
  -t {threads} \
  {params.ref} \
  {input} | \
  samtools view -hb - | \
  samtools sort \
  -@ {threads} \
  - > {output.bam} && \
  samtools index \
  -@ {threads} \
  {output.bam}
```

---

## BAMBU — ASSEMBLAGE DES TRANSCRITS

- Version : **v3.13.11**
- Conteneur au format **bambu_ebbertlab.sif** créé via Singularity à partir de l'image `library://ebbertlab/nanopore_cdna/bambu:2024-04-24` (Unique ID : `sha256.44e2b6d7282a488b95b132198b7c4ca659c9e8d6a83797493e746aa3a87ecfea`)
- Paramètres présents dans les scripts R : **bambu_guided.R** et **bambu_unguided.R** — NDR (Novel Discovery Rate) = 1 par défaut, script modifiable
- Inspirés de `bambu_guide.r` (benchmark article Nature Methods 2024) : https://github.com/WanluLiuLab/2024_LRS_AS_Benchmark_Code — https://www.nature.com/articles/s41467-024-48117-3

### Paramètres BAMBU

| Paramètre | Valeur | Description |
|---|---|---|
| `NDR` | AUTO ou 1 | Score de confiance associé à chaque transcrit assemblé |
| `min.txScore.singleExon` | Inf | Exclusion des transcrits mono-exon |
| `min.readFractionByGene` | 0.00 | Absence de filtre quantitatif |
| `trackReads` | TRUE | Génération d'un fichier reads_assignments par patient |
| `returnDistTable` | FALSE | Désactivé (coûteux en RAM) |
| `discovery` | TRUE | Découverte de nouveaux transcrits |
| `quant` | TRUE | Quantification des transcrits |

### Commande (workflow Bambu-guided / Bambu-unguided — Snakemake)

```bash
singularity exec Rscript {input.script} \
  --bam     {input.bam}        \
  --fasta   {input.genome}     \
  --output  {params.outdir}    \
  --sample  {wildcards.sample} \
  --ndr     {params.ndr}       \
  --threads {threads}
```

**Output attendu** (transcriptome par échantillon) : `{sample}.extended_annotations.gtf` — utilisé pour chaque échantillon afin de générer le catalogue complet de transcrits du Batch 0 via le module StringTie merge.

---

## FLAIR — ASSEMBLAGE DES TRANSCRITS

- Version : **v3.0.0**

### Préparation des BAM avec Intronprospector

FLAIR recommande d'utiliser Intronprospector afin d'obtenir un fichier par échantillon avec des jonctions d'épissage déduites à partir des BAM alignés (à défaut de données orthogonales de séquençage short reads couvrant les jonctions d'épissage) :

```bash
IntronProspector \
  --genome-fasta={params.ref} \
  --intron-bed6={output.bed} \
  -C 0.0 \
  {input.bam}
```

### Commande (workflow Flair_guided / Flair_unguided — Snakemake)

```bash
flair transcriptome \
  -g {params.ref} \               # génome humain
  --check_splice \
  --stringent \
  --junction_bed {input.bed} \    # fichier de nouvelles jonctions détectées — Intronprospector v1.4.0
  --junction_support 3 \          # support minimum par isoforme (défaut : 3)
  --ss_window 15 \                # fenêtre correction sites d'épissage (défaut : 15)
  --end_window 100 \              # fenêtre comparaison TSS/TES (défaut : 100)
  --max_ends 2 \                  # maximum TSS/TES par isoforme (défaut : 2)
  --no_redundant none \           # gestion des isoformes redondants (défaut : none)
  --filter default \              # filtrage des sous-isoformes (défaut : default)
  -t {threads} \
  -b {input.bam} \
  -o {params.out}
```

**Outputs attendus** (transcriptome par échantillon) :

| Fichier | Usage |
|---|---|
| `sample.isoforms.gtf` | **Seul fichier utilisé** — construction du catalogue complet du Batch 0 |
| `sample.isoforms.bed` | Coordonnées génomiques des isoformes assemblés |
| `sample.isoforms.fa` | Séquences FASTA |
| `sample.read.map.txt` | Assignements reads/transcrits |
| `sample.isoform.count.txt` | Quantification |

---

## IsoQuant — ASSEMBLAGE DES TRANSCRITS

- Version : **v3.1.0**

### Commande (workflow IsoQuant_guided / IsoQuant_unguided — Snakemake)

```bash
isoquant.py \
  --reference {input.genome} \
  --genedb {input.gtf} \        # annotation (guided mode uniquement)
  --complete_genedb \
  --bam {input.bam} \
  --data_type nanopore \
  --sqanti_output \
  --check_canonical \
  -o {params.outdir} \
  --prefix {params.prefix} \
  --threads {threads}
```

**Outputs attendus** (transcriptome par échantillon) :

| Fichier | Usage |
|---|---|
| `sample.transcript_models.gtf` | **Seul fichier utilisé** — construction du catalogue complet du Batch 0 |
| `sample.corrected_reads.bed.gz` | — |
| `sample.discovered_gene_counts.tsv` | — |
| `sample.discovered_gene_tpm.tsv` | — |
| `sample.discovered_transcript_counts.tsv` | — |
| `sample.discovered_transcript_tpm.tsv` | — |
| `sample.transcript_model_reads.tsv.gz` | — |

---

## StringTie — ASSEMBLAGE DES TRANSCRITS

- Versions : **v2.0.0** et **v3.0.3**
- Conteneur (docker image) via Singularity : `docker://quay.io/biocontainers/stringtie:3.0.3--h29c0135_0`

### Commande (workflow StringTie_guided / StringTie_unguided — Snakemake)

```bash
stringtie \
  -L \
  -p {threads} \
  -f 0.01 \
  -G {input.gtf} \    # annotation : guided mode uniquement
  -o {output} \
  {input.bam} \
  &> {log}
```

**Output attendu** (transcriptome par échantillon) : `sample.gtf` — seul fichier utilisé, servant à la construction du catalogue complet de transcrits du Batch 0.

---

## STRINGTIE MERGE — CRÉATION DU CATALOGUE COMPLET DES TRANSCRITS

### Commande (workflow RNAseq ciblé — Snakemake)

```bash
stringtie \
  -L \
  --merge \
  -p {threads} \
  -f 0 \
  -o {output.merged_gtf} \
  -G {params.transcript_guide} \
  {input.sample_gtfs}
```

---

## MODULE DE RÉALIGNEMENT DES FASTQ GUIDÉS PAR LES JONCTIONS D'ÉPISSAGE

Nouveaux BAM indexés et triés servant à la quantification :

```bash
minimap2 \
  -ax splice \
  --MD \
  -t {threads} \
  --junc-bed {input.bed} \
  {params.ref} \
  -o {output.sam} \
  {input.fastq}

samtools view \
  -b \
  -@ {threads} \
  -o {output.sorted_bam} \
  {output.sam}

samtools sort \
  -@ {threads} \
  -o {output.bam} \
  {output.sorted_bam}

samtools index \
  -@ {threads} \
  {output.bam}
```

---

## MODULE DE QUANTIFICATION / EXPRESSION — STRINGTIE

### Commande (workflow RNAseq ciblé — Snakemake)

```bash
stringtie \
  -L \
  -e \
  -p {threads} \
  -G {input.transcript_guide} \
  -o {output.output_exp} \
  -A {output.gene_abund} \
  -C {output.cov_refs} \
  {input.bam}
```

**Output** : données quantitatives — Coverage, TPM et FPKM.

---

## MODULE D'ANNOTATION — SOSTAR.py

### Commande (workflow RNAseq ciblé — Snakemake)

```bash
python3 scripts/SOSTAR.py \
  -I {params.input_folder} \
  -R {input.transcripts_ref} \
  -O {params.output_folder}
```

---

## STRINGTIE MERGE — CRÉATION DU CATALOGUE À PARTIR DES BATCHS 1 ET 2

```bash
stringtie \
  --merge \
  -L \
  -p {threads} \
  -f 0 \
  -o {output.merged_gtf} \
  {input.batch1} {input.batch2}
```

---

## MODULE DE CRÉATION DE LA MATRICE DE COMPTAGE FINALE — prepDE.py

```bash
prepDE.py \
  -i {output.sample_list} \
  -g {output.gene_counts} \
  -t {output.transcript_counts}
```

**Output** : 2 matrices de comptages créées — par gène et par transcrit.

---
---

# 🇬🇧 ENGLISH SECTION

---

## FASTQ ALIGNMENT TO THE REFERENCE GENOME (hg38) — MINIMAP2 in splice mode

```bash
minimap2 \
  -ax splice \
  --MD \
  -t {threads} \
  {params.ref} \
  -o {output.sam} \
  {input}

samtools view \
  -b \
  -@ {threads} \
  -o {output.sorted_bam} \
  {output.sam}

samtools sort \
  -@ {threads} \
  -o {output.bam} \
  {output.sorted_bam}

samtools index \
  -@ {threads} \
  {output.bam}
```

> **Note for FLAIR**: following recommendations from the official documentation; secondary alignments are excluded.

```bash
minimap2 \
  -ax splice \
  -s 80 \
  -G 200k \
  --MD \
  --secondary=no \
  -t {threads} \
  {params.ref} \
  {input} | \
  samtools view -hb - | \
  samtools sort \
  -@ {threads} \
  - > {output.bam} && \
  samtools index \
  -@ {threads} \
  {output.bam}
```

---

## BAMBU — TRANSCRIPT ASSEMBLY

- Version: **v3.13.11**
- Container format: **bambu_ebbertlab.sif**, built via Singularity from the image `library://ebbertlab/nanopore_cdna/bambu:2024-04-24` (Unique ID: `sha256.44e2b6d7282a488b95b132198b7c4ca659c9e8d6a83797493e746aa3a87ecfea`)
- Parameters defined in R scripts: **bambu_guided.R** and **bambu_unguided.R** — NDR (Novel Discovery Rate) = 1 by default (editable)
- Inspired by `bambu_guide.r` (Nature Methods 2024 benchmark): https://github.com/WanluLiuLab/2024_LRS_AS_Benchmark_Code — https://www.nature.com/articles/s41467-024-48117-3

### BAMBU Parameters

| Parameter | Value | Description |
|---|---|---|
| `NDR` | AUTO or 1 | Confidence score assigned to each assembled transcript |
| `min.txScore.singleExon` | Inf | Exclusion of single-exon transcripts |
| `min.readFractionByGene` | 0.00 | No quantitative filtering applied |
| `trackReads` | TRUE | Generates a reads_assignments file per sample |
| `returnDistTable` | FALSE | Disabled (high RAM cost) |
| `discovery` | TRUE | Novel transcript discovery enabled |
| `quant` | TRUE | Transcript quantification enabled |

### Command (Bambu-guided / Bambu-unguided workflow — Snakemake)

```bash
singularity exec Rscript {input.script} \
  --bam     {input.bam}        \
  --fasta   {input.genome}     \
  --output  {params.outdir}    \
  --sample  {wildcards.sample} \
  --ndr     {params.ndr}       \
  --threads {threads}
```

**Expected output** (per-sample transcriptome): `{sample}.extended_annotations.gtf` — used for each sample to build the complete Batch 0 transcript catalogue via the StringTie merge module.

---

## FLAIR — TRANSCRIPT ASSEMBLY

- Version: **v3.0.0**

### BAM preparation with Intronprospector

FLAIR recommends using Intronprospector to generate per-sample splice junction files derived from aligned BAMs (in the absence of orthogonal short-read sequencing data covering splice junctions):

```bash
IntronProspector \
  --genome-fasta={params.ref} \
  --intron-bed6={output.bed} \
  -C 0.0 \
  {input.bam}
```

### Command (Flair_guided / Flair_unguided workflow — Snakemake)

```bash
flair transcriptome \
  -g {params.ref} \               # human reference genome
  --check_splice \
  --stringent \
  --junction_bed {input.bed} \    # novel splice junctions from Intronprospector v1.4.0
  --junction_support 3 \          # minimum read support per isoform (default: 3)
  --ss_window 15 \                # splice site correction window (default: 15)
  --end_window 100 \              # TSS/TES comparison window (default: 100)
  --max_ends 2 \                  # maximum TSS/TES per isoform (default: 2)
  --no_redundant none \           # redundant isoform handling (default: none)
  --filter default \              # sub-isoform filtering (default: default)
  -t {threads} \
  -b {input.bam} \
  -o {params.out}
```

**Expected outputs** (per-sample transcriptome):

| File | Usage |
|---|---|
| `sample.isoforms.gtf` | **Only file used** — Batch 0 complete transcript catalogue construction |
| `sample.isoforms.bed` | Genomic coordinates of assembled isoforms |
| `sample.isoforms.fa` | FASTA sequences |
| `sample.read.map.txt` | Read-to-transcript assignments |
| `sample.isoform.count.txt` | Quantification |

---

## IsoQuant — TRANSCRIPT ASSEMBLY

- Version: **v3.1.0**

### Command (IsoQuant_guided / IsoQuant_unguided workflow — Snakemake)

```bash
isoquant.py \
  --reference {input.genome} \
  --genedb {input.gtf} \        # annotation (guided mode only)
  --complete_genedb \
  --bam {input.bam} \
  --data_type nanopore \
  --sqanti_output \
  --check_canonical \
  -o {params.outdir} \
  --prefix {params.prefix} \
  --threads {threads}
```

**Expected outputs** (per-sample transcriptome):

| File | Usage |
|---|---|
| `sample.transcript_models.gtf` | **Only file used** — Batch 0 complete transcript catalogue construction |
| `sample.corrected_reads.bed.gz` | — |
| `sample.discovered_gene_counts.tsv` | — |
| `sample.discovered_gene_tpm.tsv` | — |
| `sample.discovered_transcript_counts.tsv` | — |
| `sample.discovered_transcript_tpm.tsv` | — |
| `sample.transcript_model_reads.tsv.gz` | — |

---

## StringTie — TRANSCRIPT ASSEMBLY

- Versions: **v2.0.0** and **v3.0.3**
- Container (docker image) via Singularity: `docker://quay.io/biocontainers/stringtie:3.0.3--h29c0135_0`

### Command (StringTie_guided / StringTie_unguided workflow — Snakemake)

```bash
stringtie \
  -L \
  -p {threads} \
  -f 0.01 \
  -G {input.gtf} \    # annotation: guided mode only
  -o {output} \
  {input.bam} \
  &> {log}
```

**Expected output** (per-sample transcriptome): `sample.gtf` — only file used for Batch 0 complete transcript catalogue construction.

---

## STRINGTIE MERGE — COMPLETE TRANSCRIPT CATALOGUE CONSTRUCTION

### Command (targeted RNAseq workflow — Snakemake)

```bash
stringtie \
  -L \
  --merge \
  -p {threads} \
  -f 0 \
  -o {output.merged_gtf} \
  -G {params.transcript_guide} \
  {input.sample_gtfs}
```

---

## SPLICE JUNCTION-GUIDED FASTQ RE-ALIGNMENT MODULE

New indexed and sorted BAMs used for quantification:

```bash
minimap2 \
  -ax splice \
  --MD \
  -t {threads} \
  --junc-bed {input.bed} \
  {params.ref} \
  -o {output.sam} \
  {input.fastq}

samtools view \
  -b \
  -@ {threads} \
  -o {output.sorted_bam} \
  {output.sam}

samtools sort \
  -@ {threads} \
  -o {output.bam} \
  {output.sorted_bam}

samtools index \
  -@ {threads} \
  {output.bam}
```

---

## QUANTIFICATION / EXPRESSION MODULE — STRINGTIE

### Command (targeted RNAseq workflow — Snakemake)

```bash
stringtie \
  -L \
  -e \
  -p {threads} \
  -G {input.transcript_guide} \
  -o {output.output_exp} \
  -A {output.gene_abund} \
  -C {output.cov_refs} \
  {input.bam}
```

**Output**: quantitative data — Coverage, TPM, and FPKM.

---

## ANNOTATION MODULE — SOSTAR.py

### Command (targeted RNAseq workflow — Snakemake)

```bash
python3 scripts/SOSTAR.py \
  -I {params.input_folder} \
  -R {input.transcripts_ref} \
  -O {params.output_folder}
```

---

## STRINGTIE MERGE — CATALOGUE CONSTRUCTION FROM BATCH 1 AND BATCH 2 (GTF FILES)

```bash
stringtie \
  --merge \
  -L \
  -p {threads} \
  -f 0 \
  -o {output.merged_gtf} \
  {input.batch1} {input.batch2}
```
## 🇫🇷 MODULE D'ANNOTATION STRUCTURALE — SQANTI3 QC v5.5.3
```bash
sqanti3_qc.py \
  --isoforms {input.gtf} \
  --refGTF {input.annotation} \
  --force_id_ignore \
  --refFasta {input.genome} \
  --aligner_choice minimap2 \
  --fl_count {input.matrix_raw_counts} \
  --CAGE_peak {input.cage_peaks} \
  --polyA_peak {input.polya_sites} \
  --polyA_motif_list {input.polya_motifs} \
  --phyloP_bed {input.phyloP_bed_hboc} \
  --saturation \
  --isoform_hits \
  --report both \
  --output {params.prefix} \
  --dir {params.outdir} \
  --cpus {threads}
```
---

## FINAL COUNT MATRIX GENERATION — prepDE.py

```bash
prepDE.py \
  -i {output.sample_list} \
  -g {output.gene_counts} \
  -t {output.transcript_counts}
```
**Output**: two count matrices generated — gene-level and transcript-level.

## 🇬🇧 STRUCTURAL ANNOTATION MODULE — SQANTI3 QC Sqanti v5.5.3
```bash
sqanti3_qc.py \
  --isoforms {input.gtf} \
  --refGTF {input.annotation} \
  --force_id_ignore \
  --refFasta {input.genome} \
  --aligner_choice minimap2 \
  --fl_count {input.matrix_raw_counts} \
  --CAGE_peak {input.cage_peaks} \
  --polyA_peak {input.polya_sites} \
  --polyA_motif_list {input.polya_motifs} \
  --phyloP_bed {input.phyloP_bed_hboc} \
  --saturation \
  --isoform_hits \
  --report both \
  --output {params.prefix} \
  --dir {params.outdir} \
  --cpus {threads}
```
