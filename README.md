

# About



# Methods

## Initial inspection

$ grep -c '^>' *.fasta
cytochrome_b.fasta:108
nd1.fasta:108
prlr.fasta:12
r35.fasta:12
scar298anl.fasta:14
scar875anl.fasta:34

$ grep '^>' cytochrome_b.fasta | head -n 3
>MT795319.1 Sceloporus arenicolus isolate 2010.7 cytochrome b gene, partial cds; mitochondrial
>MT795320.1 Sceloporus arenicolus isolate CSA043 cytochrome b gene, partial cds; mitochondrial
>MT795321.1 Sceloporus arenicolus isolate CSA059 cytochrome b gene, partial cds; mitochondrial

$grep '^>' cytochrome_b.fasta | tail -n 3
>MT795424.1 Sceloporus arenicolus isolate TX2011.241 cytochrome b gene, partial cds; mitochondrial
>MT795425.1 Sceloporus arenicolus isolate TX2011.286 cytochrome b gene, partial cds; mitochondrial
>MT795426.1 Sceloporus arenicolus isolate TX2012.001 cytochrome b gene, partial cds; mitochondrial


## Getting the Outgroups

Using the outgroups from Chan et al. 2013, we get directly these following accesions

```bash
# ND1 outgroups
for acc in GQ464521 GQ464522 GQ464523 GQ464496 GQ464494 GQ464485; do
  esearch -db nucleotide -query $acc | efetch -format fasta >> outgroups_ND1_raw.fasta
done
```

```bash
# Cyt-b outgroups
for acc in GQ272800 GQ272813 AY141097 AB079242 AY141103 EU543743; do
  esearch -db nucleotide -query $acc | efetch -format fasta >> outgroups_CYTB_raw.fasta
done
```

However, the accessions from the 9 Sceloporus graciosus are not present. Neither in Chen et al. 2013, nor in Chen et al. 2020. However, we do have the voucher IDs in the figure 3 of Chet et al.2020 and the range of accessions. Using that information, I made the python script `fetch_mt_outgroups.py` to recover the proper accessions

```bash
$ python3 fetch_mt_outgroups.py
Written: ../data/sequences_by_gene/acc_outgroups_cytb.tsv
Written: ../data/sequences_by_gene/acc_outgroups_nd1.tsv

# Now, using the accesions (second column) we get the actual fasta seq

$ cut -f2 acc_outgroups_cytb.tsv | paste -sd, - \
| efetch -db nuccore -format fasta > outgroups_CYTB_sgraciosus_raw.fasta

$ cut -f2 acc_outgroups_nd1.tsv | paste -sd, - \
| efetch -db nuccore -format fasta > outgroups_ND1_sgraciosus_raw.fasta
```

Finally we check the number of seq in each file is consistent:

```bash
$ grep -c '^>' *_raw.fasta
outgroups_CYTB_raw.fasta:6
outgroups_CYTB_sgraciosus_raw.fasta:9
outgroups_ND1_raw.fasta:6
outgroups_ND1_sgraciosus_raw.fasta:9
```

A manual correction was made for:

* AB079242 from *Sceloporus occidentalis*. It actually downloaded the full mitochondrial genome, so its sequence was replaced with just the cytb gene, and the header was edited
* GQ464521 was replaced with JN6484701. Because the first was actually an error from the article, as it refered to *Urosaurus nigricaudus* instead of *Urosaurus ornatus*
* The tRNA seq was eliminated at the end of EU543743

## Tuning fasta headers

This step was performed indivdually for each fasta file using sed with regular expressions

```bash
# Example for CYTB from S. graciosus
sed 's/>.*isolate \([^ ]*\).*/>\1/' outgroups_CYTB_sgraciosus_raw.fasta > outgroups_CYTB_sgraciosus.fasta

#Example for the file with other outgroups
sed -E 's/^>[^ ]+ ([A-Z][a-z]{2})[a-z]* ([a-z]{1,6}).*/>\1_\2/' \
  outgroups_CYTB_raw.fasta > outgroups_CYTB.fasta
# Editing main from S.
sed 's/>.*isolate \([^ ]*\).*/>\1/' cytochrome_b.fasta > cytochrome_b_clean.fasta

```
## Concatenate

```bash
cat cytochrome_b_clean.fasta outgroups_CYTB.fasta outgroups_CYTB_sgraciosus.fasta > cytb_final.fasta

cat nd1_clean.fasta outgroups_ND1.fasta outgroups_ND1_sgraciosus.fasta > nd1_final.fasta
```
After joining, all dots in the header were replaced by "_"

## Aligning

```bash
# ND1
mafft --auto --reorder --adjustdirectionaccurately --thread -2 nd1_final.fasta > ../../results/phylogenetic_analysis/nd1.aln.fasta

# Cyt-b
mafft --auto --reorder --adjustdirectionaccurately --thread -2 cytb_final.fasta > ../../results/phylogenetic_analysis/cytochrome_b.aln.fasta
```
ND1: --auto chose L-INS-i (“Probably most accurate, very slow”), which uses local alignments and iteration—ideal for accuracy when you are working with a moderate number of sequences. 
123 seq × 969 pb

cyt-b: --auto chose FFT-NS-i (this is the standard strategy with iterative refinement). It is normal for this alignment size to opt for this faster option.
123 seq × 1150 pb

## Alignment trimming

For the NAD1 gene we decided to keep the alignment as it was, with very few gaps and mostly complete

```bash
# Some trials
trimal -in nd1.aln.fasta \
       -out nd1.aln.gappyout.fasta \
       -gappyout # 620 bp, seems to change too much, and that could remove signal and degrade topology (Portik et al. 2021)
trimal -in nd1.aln.fasta \
       -out nd1.aln.0.8.fasta \
       -gt 0.8 #969, doesn't change anything
```
In cytochrome b we inspected the alignment and selected the position to conserve. Our criteria was to respect the start codon and keep an alignment length multiple of 3. Thus, that left us with positions 22-1149. But in trimal pos are 0-based, so:

```bash
trimal -in cytochrome_b.aln.fasta \
       -out cytb.aln.fasta \
       -selectcols { 21-1148 } -complementary -fasta -keepheader


# Other trials
trimal -in cytochrome_b.aln.fasta \
       -out cytb.aln.gappyout.fasta \
       -gappyout    #1115 bp
trimal -in cytochrome_b.aln.fasta \
       -out cytb.aln.0.8.fasta \
       -gt 0.8    #1126 bp 
```

> ⚠️ **IMPORTANT:** MrBayes uses the same method as most maximum likelihood programs: it treats
gaps and missing characters as missing data. Thus, gaps and missing characters
will not contribute any phylogenetic information.

## Alignment concatenation

## Script cleanup.py
