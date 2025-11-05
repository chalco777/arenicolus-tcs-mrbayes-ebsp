
# I. First section
### This script reads the Genepop file, extracts the samples from each block, and assigns them to the regions as per the paper. This analysis corresponds only to microsatellite data.

# 1. Read genpop file
gen_lines <- readLines("../data/microsatellite_genotypes.gen.txt")

# 2. We find the positions where 'Pop' appears to identify the blocks
pop_idx <- grep("^\\s*[Pp][Oo][Pp]\\s*$", gen_lines)

# 3. We create a list to hold the sample IDs for each population block
pop_samples <- list()
# 4. We walk through the lines between the "Pop"s to get the IDs
for(i in seq_along(pop_idx)) {
  start_line <- pop_idx[i] + 1
  end_line   <- if(i < length(pop_idx)) pop_idx[i+1] - 1 else length(gen_lines)
  block      <- gen_lines[start_line:end_line]
  # Mantener solo las líneas con coma (cada línea de individuo tiene una coma)
  block      <- block[grepl(",", block)]
  ids <- sub("\\s*,.*$", "", block)   # eliminar todo lo que sigue a la coma
  ids <- trimws(ids)                  # eliminar espacios al inicio/final
  pop_samples[[i]] <- ids
}

# 5. We define the vector of Regions according to the paper's order
regions <- c("AA","AB","BA","BB","C","DA","DB","EA","EB","EC")

# 6. We made a df according to each sample and its region
sample_region <- data.frame(
  sample_norm = unlist(pop_samples),
  region = rep(regions, sapply(pop_samples, length)),
  stringsAsFactors = FALSE
)

# 7. Finally, we save and write the table
write.table(sample_region, file="../data/sample_to_region_microsatellite.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)



# II. Second section
## Now, we will make a tsv file that consists of a draft of sample_to_region_mtDNA.tsv, with all the samples from our mtDNA alignment file and their corresponding region (only if available from the genepop data). 

fasta_ids_raw <- sub("^>", "", grep("^>", readLines("../results/phylogenetic_analysis/ml/mtDNA_concat.fasta"), value = TRUE))
fasta_ids     <- gsub("_", "", fasta_ids_raw)

genepop_ids <- sample_region$sample


common_ids      <- intersect(fasta_ids, genepop_ids)
only_fasta       <- setdiff(fasta_ids, genepop_ids)
only_genepop     <- setdiff(genepop_ids, fasta_ids)


cat("IDs en común (FASTA ∩ Genepop):", length(common_ids), "\n")
cat("IDs solo en FASTA:", length(only_fasta), "\n")
if (length(only_fasta)) print(only_fasta)

cat("IDs solo en Genepop:", length(only_genepop), "\n")
if (length(only_genepop)) print(only_genepop)

#Base table with all fasta ids

mtDNA_all <- data.frame(
  sample      = fasta_ids_raw,
  sample_norm = fasta_ids,     # for the comparison (without _)
  stringsAsFactors = FALSE
)
print("Here")
# We join by the normalised name to bring the region
mtDNA_all <- merge(
  mtDNA_all,
  sample_region[, c("sample_norm", "region")],
  by = "sample_norm",
  all.x = TRUE
)

mtDNA_all <- mtDNA_all[, c("sample", "region")]
mtDNA_all <- mtDNA_all[order(mtDNA_all$region), ]

write.table(mtDNA_all, file = "../data/sample_to_region_mtDNA.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Note that the rest of the samples without region where filled using the Figure 3 from Chan et al. 2020. The count of each sample was added based on that figure too.
# The final manually edited file is  ""../data/sample_to_region_mtDNA.tsv", with the same name