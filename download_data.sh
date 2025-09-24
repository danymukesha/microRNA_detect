#!/usr/bin/env bash
set -euo pipefail

mkdir -p data
cd data

echo "1) Download chr19 VCF (1000 Genomes phase 3)"
VCF_URL="https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
wget -c "$VCF_URL" -O ALL.chr19.phase3.vcf.gz

echo "2) (optional) download index .tbi if present"
wget -c "${VCF_URL}.tbi" -O ALL.chr19.phase3.vcf.gz.tbi || true

echo "3) Download chr19 FASTA from UCSC"
# UCSC chromosomes directory contains chr19.fa.gz; i used chromFa archive, but i can also use per-chrom
CHR19_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz" # I used hg19 here because 1000 Genomes phase 3 is natively listed in hg19 coordinates at the UCSC mirror
wget -c "$CHR19_URL" -O chr19.fa.gz
gunzip -f chr19.fa.gz

echo "4) Download GENCODE v19 GTF (annotation for hg19)"
#GTF_URL="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
wget -c "$GTF_URL" -O gencode.v19.annotation.gtf.gz
gunzip -f gencode.v19.annotation.gtf.gz

echo "5) Download miRBase mature sequences"
#MATURE_URL="https://mirbase.org/ftp/CURRENT/mature.fa.gz"
MATURE_URL="https://www.mirbase.org/download/mature.fa"
wget -c "$MATURE_URL" -O mature.fa.gz
gunzip -f mature.fa.gz

echo "All downloads complete. Files in data/:"
ls -lh
