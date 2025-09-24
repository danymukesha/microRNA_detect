# Method

> **Methods (short)**
> We developed a lightweight, reproducible pipeline to identify candidate regulatory variants that alter predicted microRNA (miRNA) seed matches in human genes of interest. We downloaded 1000 Genomes Phase 3 variant calls (per-chromosome VCFs) and extracted 3′-UTR coordinates from GENCODE v19 (GRCh37/hg19). For each variant overlapping a 3′-UTR we extracted a ±30 bp local sequence from the matching reference genome and constructed the alternate allele sequence. We compiled mature human miRNA sequences from miRBase and derived canonical 7-nt seed sequences (positions 2–8). For each variant we tested presence/absence of the reverse-complement seed motif in the reference and alternate sequences. Variants that eliminate a seed were annotated as predicted **seed loss**; variants that create a seed were annotated as predicted **seed gain**. We recorded the number and identity of miRNAs gained and lost and computed an `impact` score equal to the sum of gains and losses. We prioritized variants by impact, allele frequency, and concordance with external resources. Candidate variants were cross-referenced with dbSNP/ClinVar (NCBI E-utilities) and Ensembl VEP for consequence annotation; colocalization with expression quantitative trait loci (eQTLs) was checked using the GTEx portal API. All scripts used in the pipeline are available at \[repository link], and the analysis can be run on a personal computer with Python (cyvcf2, Biopython, pandas). The code and input URLs are provided in the Supplementary Materials.

# Quick top-candidate summary (for APOE and TOMM40)

During my analysis, i generated a file and ranked candidates by `impact = n_GAIN + n_LOSS`. Here are the **top 5** (columns simplified). Run the annotation script below on these rsIDs to get dbSNP/ClinVar/GTEx status.

1. CHR: 19 — POS: 45,410,250 — REF>ALT: G>C — ID: **rs558171988** — GENE: APOE — n\_GAIN=14 — n\_LOSS=7 — impact=21
2. CHR: 19 — POS: 45,411,748 — REF>ALT: C>T — ID: **rs568735371** — GENE: APOE — n\_GAIN=12 — n\_LOSS=9 — impact=21
3. CHR: 19 — POS: 45,402,516 — REF>ALT: G>A — ID: **rs113112231** — GENE: APOE — n\_GAIN=5 — n\_LOSS=12 — impact=17
4. CHR: 19 — POS: 45,409,089 — REF>ALT: C>T — ID: **rs553197505** — GENE: APOE — n\_GAIN=5 — n\_LOSS=12 — impact=17
5. CHR: 19 — POS: 45,394,544 — REF>ALT: C>T — ID: **rs554906988** — GENE: TOMM40 — n\_GAIN=11 — n\_LOSS=5 — impact=16

(the full top-20 table will be saved later to `/data/top20_candidates.csv` on the gh repo of this tool.)

These variants are predicted to gain or lose multiple human miRNA seeds in their local 61-bp windows and therefore are prioritized. 
Some have rsIDs already (which makes cross-referencing straightforward).

# B Variant prioritization ranking template

Columns (header names + description):

* `CHROM` — chromosome (e.g., 19)
* `POS` — 1-based genomic coordinate
* `REF` — reference allele
* `ALT` — alternate allele
* `ID` — rsID if present (e.g., rs12345)
* `GENE` — overlapping gene(s) from GTF (e.g., APOE)
* `FEATURE` — region type (3'UTR / 5'UTR / exon / intron / intergenic)
* `STRAND` — gene strand (+/-)
* `AF_1000G` — allele frequency from 1000 Genomes (if available)
* `AF_gnomAD` — allele frequency from gnomAD (optional add-on via API)
* `n_GAIN` — number of miRNAs predicted to **gain** binding (from seed RC matches)
* `n_LOSS` — number of miRNAs predicted to **lose** binding
* `impact` — `n_GAIN + n_LOSS` (primary ranking metric)
* `seed_types` — list or pipe-delimited list of gained/lost miRNA names (e.g., GAIN\:miR-..|LOSS\:miR-..)
* `consequence` — VEP/Ensembl predicted consequence (e.g., 3\_prime\_UTR\_variant)
* `ClinVar_status` — ClinVar aggregate classification (Pathogenic / Benign / VUS / Not found)
* `dbSNP_exists` — yes/no and dbSNP rsID (if different)
* `GTEx_eQTL` — yes/no (+ tissue(s) if present) — points to colocalization with expression QTLs for the gene
* `cross_species_conserved` — PhyloP/PhastCons score (optional; indicates conservation)
* `prioritization_score` — combined score (custom; see suggested formula below)
* `notes` — free text for manual curation

Suggested simple prioritization score (0–100):
`prioritization_score = clamp( 50*(impact/max_impact) + 30*(1 - log10(AF+1e-6)/log10(0.5)) + 20*(1 if ClinVar_pathogenic else 0) , 0,100)`

> High impact gets weight; very rare variants get a boost (rare = possibly higher effect), ClinVar pathogenicity heavily weights up. I tried to tweak weights for my objectives (e.g., prioritize common regulatory variants vs. rare variants).


# Figure legend (for a single-panel figure that visualizes prioritization) [#I have to upload it later#]

Figure: miRNA-seed disruption prioritization for variants in the APOE locus.

Legend text:

> **Figure 2. Candidate regulatory variants in the APOE/TOMM40 locus predicted to alter miRNA binding.** (a) Manhattan-like plot across the APOE locus showing the `impact` score (n\_gain + n\_loss) per variant (y-axis) against genomic position (x-axis). Variants are colored by allele frequency (heatmap: common → rare). Top candidates (labeled by rsID) are highlighted and annotated with the number of miRNAs gained and lost. (b) Example sequence-level depiction for the top candidate (rs558171988) showing the reference 61-nt window, alternate sequence, and the gained/lost miRNA seed matches (7-mer) in colored boxes. (c) Prioritization table excerpt showing variant, gene, impact, AF (1000G), ClinVar status, and GTEx eQTL overlap. The pipeline integrates population frequency and public annotations to rank candidates for follow-up functional assays (e.g., 3′-UTR luciferase with miRNA mimic).


