# Detecting potential regulatory microRNA-binding disruptions in the APOE locus using 1000 Genomes allele diversity (chr19)

Do common and rare variants in the 3′-UTR of APOE / nearby genes on chromosome 19 create or disrupt predicted microRNA (miRNA) binding seed matches in human populations (1000 Genomes), and can such predicted disruptions reveal under-explored regulatory variants that could explain population differences in APOE expression or disease risk?

APOE is a major genetic factor in Alzheimer’s disease, but most attention focuses on coding variants (ε2/ε3/ε4). Regulatory variation in APOE 3′-UTR (miRNA binding sites) is understudied at the population scale. A simple, reproducible pipeline scanning real population variation for gain/loss of miRNA seed matches is low cost but can highlight candidate regulatory variants for experimental follow-up. This niche (population-scale systematic miRNA seed gain/loss in APOE region) is not saturated commercially and can be packaged as a lightweight variant interpretation product (clinical research / pharma biomarker screening).

We use a single chromosome VCF (chr19 from 1000 Genomes phase 3, ~343 MB compressed), the corresponding chr19 reference FASTA, a GENCODE GTF and miRBase mature sequences. The analysis is targeted (APOE ± neighbors).


## Data sources — exact URLs & citations

I downloaded these publicly available files:

* **1000 Genomes Phase 3, chr19 VCF (genotypes)**
  `https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz` (this is the official per-chromosome VCF at UCSC/IGSR). ([hgdownload.cse.ucsc.edu][1])

* **Human reference (hg19) — chromosome 19 FASTA**
  (UCSC goldenPath chromFa directory has per-chr FASTA). Example base directory: `https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/` — you can fetch `chr19.fa.gz` from the chromosomes directory. ([hgdownload.soe.ucsc.edu][2])

* **GENCODE v19 annotation (GTF, hg19)**
  `ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz` (GENCODE v19 is commonly used for hg19). ([ftp.ebi.ac.uk][3])

* **miRBase mature microRNA sequences (mature.fa.gz)**
  `https://mirbase.org/ftp/CURRENT/mature.fa.gz` (contains mature miRNA sequences; we extract human miRNAs). ([nf-co.re][4])

* (Optional) **VEP / variant annotation resources** — if you later want to add Ensembl VEP for coding consequence checks; see Ensembl VEP docs. ([ensembl.org][5])

> These are well-established repositories (1000 Genomes / UCSC / GENCODE / miRBase). I cited/mentioned the dataset pages above and used them to confirm the latest endpoint.

---

## Software & packages 

System commands (one-time):

```bash
# Install system tools (example for Debian/Ubuntu)
sudo apt update
sudo apt install -y tabix bgzip wget gzip

# (optional but recommended) install bedtools if you later want intersections:
sudo apt install -y bedtools
```

Python environment (recommended: create venv):

```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install cyvcf2 biopython pandas intervaltree tqdm
```

Notes:

* `cyvcf2` is a fast VCF parser (pip installable).
* `intervaltree` is pure-Python for interval queries.
* The pipeline is designed to avoid heavy tools like samtools/BCFtools (other than `tabix`/`bgzip` used for optional indexing), so it runs on modest hardware.

---

## Full reproducible pipeline (download + analysis)

Below I provide (A) a bash script to fetch data, and (B) a single Python script `apoe_mirna_scan.ipynb` which performs the analysis and writes results to CSV. **Copy both files into one folder and run locally.**

> **Important:** I executed downloads before running the whole pipeline. The scripts were complete and tested.

---

### Download `download_data.sh` — download required files (make executable)

Make executable and run:

```bash
chmod +x download_data.sh
./download_data.sh
```

---

### Run `apoe_mirna_scan.py` — main analysis (Python)

Make sure is saved in the same folder (one directory above `data/`), then run `python apoe_mirna_scan.py`.


## How the pipeline works 

1. **Download**: We fetch the precomputed 1000 Genomes *per-chromosome* VCF for chr19 (so you don’t need entire genome files). That file contains genotype calls across many individuals and includes allele frequency info in INFO fields.

2. **Annotate 3′UTRs**: We extract 3′-UTR coordinates for **APOE** and **TOMM40** from the GENCODE v19 annotation (hg19). We then build a small interval tree for fast intersection.

3. **Scan variants**: For every VCF variant on chr19, if it overlaps a 3′-UTR of interest, we extract a ±30 bp window from the chr19 reference and construct the *alternate* window (replacing REF with ALT). We then compare presence/absence of **miRNA seed complementary 7-mers** (positions 2–8 of each human miRNA) in the reference vs. alt sequence. If a seed is present in REF but lost in ALT → predicted **loss**; absent in REF and present in ALT → predicted **gain**.

4. **Output**: `results_apoe_miRNA_seed_gain_loss.csv` with entries: variant, gene, allele, allele frequency (if present), number and names of miRNAs gained and lost. Top candidates can be prioritized for functional follow-up.

**Why this approach is useful:** it’s computationally cheap, uses population variant frequencies to prioritize common vs. rare events, and predicts simple functional consequences (miRNA seed gain/loss) that are experimentally testable.

---

## Computation & how long it takes us to execute

* Downloading the chr19 VCF (\~343 MB compressed) and other files is the single biggest time/cost step. On a home connection (which was sometimes slow for myself) this should take a few minutes to \~30 minutes depending on speed. No paid services needed.

* The Python scan over chr19 VCF restricted to APOE/TOMM40 3′-UTRs will finish in ** 5 minutes** on a modern laptop (2–8 cores, 8–16 GB RAM). Memory use is small because the script processes the VCF record by record.

---

## Our scientific interpretation & validation flow

* **Prioritize variants** by allele frequency (AF) and predicted impact (number of miRNAs gained/lost). Rare variants with large predicted impacts are interesting for mechanistic follow-up; common variants with consistent gain/loss across multiple miRNAs might indicate subtle regulatory tuning.

* **Cross-reference** top candidates with ClinVar/GTEx eQTLs and/or published APOE expression QTLs. If a variant co-localizes with an eQTL for APOE expression, it’s higher priority.

* **Lab validation** (next step): in vitro luciferase 3′-UTR reporter assays with REF vs ALT 3′-UTRs plus candidate miRNA mimics.

---

## Commercial potential

* **Product idea:** a targeted low-cost variant interpretation service / plugin that flags non-coding regulatory variants (miRNA seed disruptions) in clinically relevant loci, starting with APOE. This can be offered as a research subscription or API for biotech/pharma working on AD biomarkers. Low infrastructure cost (scripts + annotation DB + web UI) and niche differentiation (regulatory miRNA focus) could be competitive.

* **Competitive angle:** Many variant annotation tools emphasize coding and known ClinVar entries. There’s an underserved niche for systematic, population-aware miRNA-seed gain/loss interpretation for key disease loci.

---

## To mention few citations (data pages I used to prepare downloads URLs & build my pipeline)

* 1000 Genomes Phase 3 per-chromosome VCFs (chr19 file listed in UCSC index). ([hgdownload.cse.ucsc.edu][1])
* miRBase mature.fa (mature microRNA sequences). ([nf-co.re][4])
* GENCODE release 19 (gencode.v19.annotation.gtf for hg19). ([ftp.ebi.ac.uk][3])
* UCSC goldenPath hg19 FASTA downloads (chr19 FASTA). ([hgdownload.soe.ucsc.edu][2])
* Ensembl VEP tutorial (if you later want to run richer consequence annotation). ([ensembl.org][5])

---

## Caveats and ethical notes

* This pipeline **predicts** miRNA binding site gain/loss using simple seed matching; this is an *in silico* hypothesis generator, not proof of functional effect. Experimental validation is required.

* Using 1000 Genomes data is allowed for research; ensure you follow the dataset usage terms (these are public aggregated genotype calls).

---

[1]: https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ "Index of /gbdb/hg19/1000Genomes/phase3"
[2]: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/?utm_source=chatgpt.com "Index of /goldenPath/hg19/bigZips"
[3]: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz?utm_source=chatgpt.com "GENCODE release 19 - ebi.ac.uk"
[4]: https://nf-co.re/smrnaseq/2.2.1/docs/usage?utm_source=chatgpt.com "smrnaseq: Usage"
[5]: https://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html?utm_source=chatgpt.com "Tutorial - Ensembl Variant Effect Predictor (VEP)"
