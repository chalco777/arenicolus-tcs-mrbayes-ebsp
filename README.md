

# About



# Methods

## Initial inspection

```bash
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
```

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

## Alignment trimming and concatenation

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

Finally, we concatenate the alignments using AMAS. With the following command we get the final alignment in fasta and nexus, as well as the necessary partitions file for running ModelFinder in iqtree

```bash
python -m amas.AMAS concat \
  -i cytb.aln.fasta nd1.aln.fasta \
  -f fasta -d dna \
  -u fasta -t mtDNA_concat.fasta \ #ran a second time for the second output
  -u nexus -t mtDNA_concat.nex \
  -p mtDNA_partitions_by_gene.txt
```
> **Note:** You can specify the output partition format with `--part-format raxml`

## Model selection

As done in Chan et al. 2020, we define 6 diferent partitions, by codon position and by gene.

```bash
echo 'DNA, cytb_pos1 = 1-1128\3
DNA, cytb_pos2 = 2-1128\3
DNA, cytb_pos3 = 3-1128\3
DNA, nd1_pos1  = 1129-2097\3
DNA, nd1_pos2  = 1130-2097\3
DNA, nd1_pos3  = 1131-2097\3' > mtDNA_partitions_codon.txt
```

Next, we use iqtree's implementation of ModelFinder to select the best models for each partition, using BIC as our criterion

```bash
iqtree2 --prefix mtDNA \
        -s mtDNA_concat.fasta \
        -p mtDNA_partitions_codon.txt \
        -m MFP --merit BIC -T AUTO \
        -mset mrbayes

# Another tested

iqtree2 --prefix mtDNA \
        -s mtDNA_concat.fasta \
        -p mtDNA_partitions_codon.txt \
        -m MFP --merit AIC -T AUTO \
        -mset mrbayes

```
* -m for ModelFinder and tree inference
* BIC penalizes more extra params than AIC
* mset mrbayes restricts models to those present in MrBayes

The selected best model for each partition, together with the commands for Mr Bayes, were added at the end of [mtDNA_concat.nex](results/phylogenetic_analysis/model_selection/mtDNA_concat.nex)

```bash
mb
```

> **Note:** In the context of explaining the analysis of primates.nexus mitochondrial seq, MrBayes manual mentions: The Code setting is only relevant if the Nucmodel is set to Codon. The Ploidy setting is also irrelevant for us. 

## MrBayes and ML


Using ML

```bash
iqtree2 \
  --prefix mtDNA_ML \
  -s mtDNA_concat.fasta \
  -p mtDNA_partitions_codon.txt \
  -m MFP --merit BIC \
  -B 1000 --bnni \
  -T AUTO -seed 20251104

```
Revisar log

## DNAsp and R diversity stats:

For these analyses only individuals present in the GenPop file (that is, those with microsatellite data) had the specific sample region. The regions for each sample were extracted from the microsatellite_genotypes.gen.txt file using the (first section of the) script [sample_to_regions.R](scripts/sample_to_regions.R). For individuals who do not appear in Genepop, there is no public list in the supplementary materials that assigns their regions. These are likely to be:

* Individuals sequenced only for mtDNA (without microsatellite genotypes), including several outgroups such as Phr_corona, museum specimens (CAS, MVZ, MSB, etc.), and old catches.
* Individuals with no location information or incomplete data, which were therefore excluded from the structure analyses.

It is worth noting that for calculating S, pi, number of haplotypes and K, from the concatenated mitochondrial alignment they used the previous, not deduplicated alignment. Thus, to try to replicate the original DNAsp stats from the concatenated alignment, we first generated the [sample_to_region_mtDNA.tsv](data/sample_to_region_mtDNA.tsv) file using the second section of [sample_to_regions.R](scripts/sample_to_regions.R). This file only had the regions for the sample ids that also happened to be sequenced for microsatellites. Then we manually annotated the count for each sample id (or voucher) from the Figure 3 in Chan et al. 2020. We also tried to complete the regions that were not given using that phylogenetic tree. 

After that, we used the count for each mitochondrial sequence (or haplotype) together with the script [reverse_collapse_cleanup.py](scripts/reverse_collapse_cleanup.py) for regenerated their original, not deduplicated, alignment. 

```bash
python reverse_collapse_cleanup.py \
    -i ../results/phylogenetic_analysis/alignment/mtDNA_concat.nex \
    --counts ../data/sample_to_region_mtDNA.tsv \
    --output ../results/haplotypes/diversity_stats/mtDNA_concat_reversed.nex
```
Our new matrix has 223x2097 bp, in contrast with theirs, that had 225x2097 bp. Hower, it's much more close than our previous, deduplicated matrix of 123x2097 bp. We will compare both, our dedup and our regenerated original, with the original, not deduplicated, matrix of Chan.

For diversity stats, we removed outgroups with AMAS:

```bash
TAXA="JWA338 JWA470 Phr_corona SGR4 SGR3 Sce_jarrov Sce_merria Sce_occide Uro_ornatu Uta_stansb MVZ149956 MVZ237413 MVZ241596 CAS223822 CAS229140"

amas remove -d dna -f nexus -i mtDNA_concat.nex -x $TAXA -u nexus -g tmp_
mv tmp_mtDNA_concat.nex-out.nex alignment_filtered/mtDNA_concat_filtered.nex

# mtDNA_concat_reversed.nex -> NEXUS filtrado
amas remove -d dna -f nexus -i mtDNA_concat_reversed.nex -x $TAXA -u nexus -g tmp_
mv tmp_mtDNA_concat_reversed.nex-out.nex alignment_filtered/mtDNA_concat_reversed_filtered.nex
```

We used DnaSP and the R packages ape and pegas for calculating diversity statistics. The R notebook used is [diversity_stats.Rmd](/scripts/diversity_stats.Rmd). As the default algorithm of DnaSP does, in our R script we ignored columns with gaps and missing nucleotides before calculating the metrics. Preliminar results (yet to be commented) are [here](/results/haplotypes/diversity_stats/).

## PopArt

We prepared the input nexus with the traits section using the python script [popart_region_nexus.py](/scripts/popart_region_nexus.py)

```bash

python popart_region_nexus.py ../data/sample_to_region_mtDNA.tsv ../results/haplotypes/diversity_stats/alignment_filtered/mtDNA_concat_reversed_filtered.nex ../results/haplotypes/popart/mtDNA_concat_reversed_filtered_traits.nex
```


According to PopArt: # of parsimony informative sites 65


> **Note:** PopART masks (ignores) any alignment column containing gaps or ambiguous characters (?, N, Y, R) before collapsing. This means that columns with any gaps/ambiguities are removed for everyone, and some sequences that were not identical may become identical after this masking (and therefore collapse). In this sense PopArt behaves similarly to DnaSP's default mode

(Luego replicar para loci nucleares)

Count number of PI sites. They used PAUP. I will use iqtree result

mtDNA_concat.fasta->They count # of Haplotypes, S (segregating sites), pi (nucleotide div) and k (avg. pairwise diff). I will use DivGenSeq.R (or ape+pegas)

Redes de haplotipos TCS network
Alineamiento->Debo generar un nex con haplotipos y con los tags y de donde viene aca haplotipo ps, o sea a qué tag pertence. Parece que ella lo genero desde DNAsp ese nex
