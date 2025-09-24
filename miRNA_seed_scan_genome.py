#!/usr/bin/env python3
"""
Scans variants in a VCF for predicted miRNA seed gain/loss events across a
user-specified genomic region or whole genome. Supports GRCh37 (hg19) and GRCh38.

Usage:
    python miRNA_seed_scan_genome.py \
        --vcf data/ALL.chr19.phase3.vcf.gz \
        --ref data/chr19.fa \
        --gtf data/gencode.v19.annotation.gtf \
        --mature data/mature.fa \
        --assembly hg19 \
        --genes APOE,TOMM40 \
        --out results.csv

For whole-genome scanning, provide a genome VCF (or loop over per-chrom VCFs).

Notes & usage:
For whole-genome scanning: provide a whole-genome VCF (e.g., ALL.chr1-22.phase3... or 
loop through per-chromosome files). If scanning all chromosomes, ensure your --ref FASTA 
includes all chromosomes in the same naming convention as your GTF (e.g., chr1, chr2, ...
vs 1,2,...). GENCODE GTF v19 is for hg19; for GRCh38 use matching GENCODE release and FASTA.

For GRCh38: download GENCODE v29+ or appropriate GRCh38 GTF and the GRCh38 chromosome FASTA;
pass --assembly GRCh38 (script not changing behavior by assembly but you must supply matching inputs).

Execution time: whole-genome will take proportionally longer; streaming the VCF keeps memory low
"""
import argparse
import os
import sys
from cyvcf2 import VCF
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from intervaltree import IntervalTree
from tqdm import tqdm

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

def load_miRNA_seeds(mature_fa_path):
    seeds = {}
    for rec in SeqIO.parse(mature_fa_path, "fasta"):
        header = rec.id
        if header.startswith("hsa-") or header.startswith("hsa"):
            seq = str(rec.seq).upper()
            if len(seq) >= 8:
                seed = seq[1:8]
                seeds[header] = seed
    return seeds

def parse_gtf_three_prime_utrs(gtf_path, genes_of_interest=None):
    gene_utrs = {}
    with open(gtf_path,'rt') as fh:
        for line in fh:
            if line.startswith("#"): continue
            parts = line.split("\t")
            #if len(parts) < 9: continue
            chrom, src, feature, start, end, score, strand, frame, attr = parts
            if feature != "three_prime_utr": continue
            # parse gene_name
            attrs = {}
            for kv in attr.strip().split(";"):
                kv = kv.strip()
                if kv=="":
                    continue
                if " " in kv:
                    k,v = kv.split(" ",1)
                    attrs[k] = v.strip().strip('"')
            gene_name = attrs.get("gene_name") or attrs.get("gene_id")
            if (genes_of_interest is None) or (gene_name in genes_of_interest):
                gene_utrs.setdefault(chrom, []).append((int(start), int(end), gene_name, strand))
    return gene_utrs

def build_interval_trees(utr_dict):
    trees = {}
    for chrom, segs in utr_dict.items():
        it = IntervalTree()
        for (s,e,g,strand) in segs:
            it[s-1:e] = (g,strand)
        trees[chrom] = it
    return trees

def seq_window(chrom_seq, pos1, window):
    L = len(chrom_seq)
    start = max(1, pos1 - window)
    end = min(L, pos1 + window)
    return chrom_seq[start-1:end]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Input VCF (can be whole genome or chr-specific)")
    ap.add_argument("--ref", required=True, help="Reference FASTA for matching assembly (per-chrom or genome)")
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (matching assembly)")
    ap.add_argument("--mature", required=True, help="miRBase mature.fa")
    ap.add_argument("--assembly", choices=("hg19","GRCh37","hg38","GRCh38"), default="hg19")
    ap.add_argument("--genes", default="", help="Comma-separated gene list to restrict to (optional)")
    ap.add_argument("--window", type=int, default=30, help="Window size to search around variant (bp)")
    ap.add_argument("--out", default="results_genome_miRNA.csv")
    args = ap.parse_args()

    genes = [g.strip() for g in args.genes.split(",") if g.strip()] if args.genes else None

    print("Loading miRNA seeds...")
    seeds = load_miRNA_seeds(args.mature)
    seed_rc_map = {}
    for name, seed in seeds.items():
        rc = revcomp(seed)
        seed_rc_map.setdefault(rc, []).append(name)

    print("Parsing GTF for 3'UTR segments...")
    utrs = parse_gtf_three_prime_utrs(args.gtf, genes_of_interest=genes)
    trees = build_interval_trees(utrs)

    # load reference FASTA; could be per-chrom file, but here we expect a single fasta with chrom id matching VCF chroms
    ref_seqs = {}
    for rec in SeqIO.parse(args.ref, "fasta"):
        rid = rec.id
        ref_seqs[rid] = str(rec.seq).upper()

    print("Streaming VCF:", args.vcf)
    vcf = VCF(args.vcf)
    results = []
    for var in tqdm(vcf):
        chrom = str(var.CHROM)
        pos = var.POS
        # normalize chrom keys to match FASTA & GTF ids: try variants
        tree = None
        for key in (chrom, "chr"+chrom, chrom.replace("chr","")):
            if key in trees:
                tree = trees[key]
                chrom_key = key
                break
        if tree is None:
            # no UTRs in this chrom (skip)
            continue
        hits = tree[pos-1]
        if not hits:
            continue
        ref = var.REF
        alts = var.ALT
        af = var.INFO.get('AF') if 'AF' in var.INFO else None
        # get chrom_seq
        chrom_seq = None
        for key in (chrom, "chr"+chrom):
            if key in ref_seqs:
                chrom_seq = ref_seqs[key]
                break
        if chrom_seq is None:
            # skip if we don't have reference sequence for this chrom
            continue

        for alt in alts:
            seq_ref_win = seq_window(chrom_seq, pos, args.window)
            left_len = min(args.window, pos-1)
            center_idx = left_len
            left = seq_ref_win[:center_idx]
            right = seq_ref_win[center_idx+len(ref):]
            alt_seq_win = left + alt + right

            gained = []
            lost = []
            for rc_seed, mirnas in seed_rc_map.items():
                if rc_seed in seq_ref_win and rc_seed not in alt_seq_win:
                    lost.extend(mirnas)
                if rc_seed not in seq_ref_win and rc_seed in alt_seq_win:
                    gained.extend(mirnas)
            gained_set = sorted(set(gained))
            lost_set = sorted(set(lost))
            results.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "GENE": ";".join(sorted(set([h.data[0] for h in hits]))),
                "STRAND": ";".join(sorted(set([h.data[1] for h in hits]))),
                "AF": af,
                "GAIN_miRNAs": "|".join(gained_set),
                "LOSS_miRNAs": "|".join(lost_set),
                "n_GAIN": len(gained_set),
                "n_LOSS": len(lost_set),
                "ID": var.ID if var.ID else "."
            })

    df = pd.DataFrame(results)
    if not df.empty:
        df['impact'] = df['n_GAIN'] + df['n_LOSS']
    df.to_csv(args.out, index=False)
    print("Wrote", args.out)

if __name__ == "__main__":
    main()
