library(ape)
library(ggtree)
install.packages("ggtree")
# Instalar BiocManager si no lo tienes
install.packages("BiocManager")
install.packages('gdtools')
install.packages('ggiraph')
# Cargar la librería
library(ggtree)

# We read the final ML tree
tree <- read.tree("results/phylogenetic_analysis/ml/mtDNA_ML.treefile")

outgroups <- c("JWA338", "JWA470", "Phr_corona", "SGR4", "SGR3", 
               "Sce_jarrov", "Sce_merria", "Sce_occide", "Uro_ornatu", 
               "Uta_stansb", "MVZ149956", "MVZ237413", "MVZ241596", 
               "CAS223822", "CAS229140")


