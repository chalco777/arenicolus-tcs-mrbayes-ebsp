# ============================================================
# ANALISIS DE STRUCTURE: QUÉ POBLACIONES CORRESPONDEN A CADA CLUSTER
# ============================================================

library(tidyverse)
library(stringr)
library(viridis)
library(here)

# Asume que se abre el proyecto (.Rproj) desde RStudio, por lo que getwd() ya es la raíz
results_dir <- file.path("results", "structure")
meta_file <- file.path("data", "fst_analysis", "lagartijas_STRUCTURE.txt")

if (!dir.exists(results_dir)) {
  stop("No existe el directorio de resultados generado por run_structure_parallel.sh: ", results_dir)
}

if (!file.exists(meta_file)) {
  stop("No se encuentra el archivo de metadatos de STRUCTURE: ", meta_file)
}

# ============================================================
# 1. LEER METADATA DESDE ARCHIVO ORIGINAL DE STRUCTURE
# ============================================================

meta <- read.table(meta_file,
                   header = FALSE, fill = TRUE, stringsAsFactors = FALSE)

colnames(meta)[1:2] <- c("Ind", "POP")
meta <- meta %>% select(Ind, POP)

cat("📁 Metadata cargada:", nrow(meta), "individuos.\n")

# ============================================================
# 2. LEER LOG-PROBABILIDADES (LnP) DE LOS ARCHIVOS _f
# ============================================================

files <- list.files(results_dir, pattern = "_f$", full.names = TRUE)
if (length(files) == 0) {
  stop("No se encontraron archivos _f en: ", results_dir)
}

extract_lnp <- function(path) {
  txt <- readLines(path, warn = FALSE)
  lnL_line <- txt[grep("Estimated Ln Prob of Data", txt)]
  lnL_val <- as.numeric(str_extract(lnL_line, "-?\\d+\\.\\d+"))
  
  fname <- basename(path)
  K <- as.numeric(str_extract(fname, "(?<=K)\\d+(?=_)"))
  rep <- as.numeric(str_extract(fname, "(?<=rep)\\d+"))
  
  tibble(K = K, Replicate = rep, LnP = lnL_val)
}

df <- map_dfr(files, extract_lnp)

cat("📊 Runs cargados:", nrow(df), "\n")

# ============================================================
# 3. CALCULAR EL MEJOR K (MÉTODO DE EVANNO)
# ============================================================

evanno <- df %>%
  group_by(K) %>%
  summarise(
    Mean_LnP = mean(LnP),
    SD_LnP = sd(LnP),
    .groups = "drop"
  ) %>%
  arrange(K)

evanno$LnP_prime <- c(NA, diff(evanno$Mean_LnP))
evanno$LnP_doubleprime <- c(NA, diff(evanno$LnP_prime[-1]), NA)
evanno$DeltaK <- abs(evanno$LnP_doubleprime) / evanno$SD_LnP

bestK <- evanno %>% filter(DeltaK == max(DeltaK, na.rm = TRUE)) %>% pull(K)

cat("✨ Mejor K según Evanno:", bestK, "\n")

# ============================================================
# 4. FUNCIÓN PARA EXTRAER Q-MATRIX DE UN ARCHIVO _f
# ============================================================

extract_Qmatrix <- function(path) {
  txt <- readLines(path, warn = FALSE)
  
  # Encontrar inicio de Q-matrix
  start <- grep("Inferred ancestry of individuals", txt)
  if (length(start) == 0) return(NULL)
  start <- start + 2
  
  # Fin: primera línea vacía luego de la lista
  rel <- which(txt[(start+1):length(txt)] == "")[1]
  end <- if (is.na(rel)) length(txt) else start + rel - 1
  
  block <- txt[start:(end-1)]
  
  # Leer tabla cruda
  df_raw <- read.table(text = block, fill = TRUE, stringsAsFactors = FALSE)
  
  # LIMPIAR ASTERISCOS
  df_raw[] <- lapply(df_raw, function(x) gsub("\\*", "", x))
  
  # CONVERTIR A NUMÉRICOS DONDE APLICA
  suppressWarnings(df_raw[] <- lapply(df_raw, function(x) type.convert(x, as.is = TRUE)))
  
  # Ahora df_raw tiene columnas así:
  # Index | Ind | POP? | Q1 | Q2 | ... | QK
  
  # Identificar columnas Q verdaderas (probabilidades entre 0 y 1)
  is_prob <- function(col) is.numeric(col) && all(col >= 0 & col <= 1, na.rm = TRUE)
  
  prob_cols <- sapply(df_raw, is_prob)
  
  df <- df_raw[, c("V1", "V2", names(df_raw)[prob_cols])]
  
  colnames(df)[1:2] <- c("Index", "Ind")
  
  # renombrar Q
  K <- sum(prob_cols)
  colnames(df)[3:(2+K)] <- paste0("Q", 1:K)
  
  return(df)
}


# ============================================================
# 5. SELECCIONAR LA MEJOR RÉPLICA DEL MEJOR K
# ============================================================

files_bestK <- list.files(results_dir, pattern = paste0("K", bestK, "_rep.*_f$"), full.names = TRUE)

best_rep <- df %>%
  filter(K == bestK) %>%
  arrange(desc(LnP)) %>%
  slice(1) %>%
  pull(Replicate)

best_file <- files_bestK[grep(paste0("rep", best_rep), files_bestK)]
if (length(best_file) == 0) {
  stop("No se encontró el archivo _f para K=", bestK, " réplica ", best_rep, " en ", results_dir)
}
best_file <- best_file[1]
Qmat <- extract_Qmatrix(best_file)

cat("🧪 Usando réplica óptima:", best_rep, "para K =", bestK, "\n")

# ============================================================
# 6. UNIR Q-MATRIX + POP
# ============================================================

Qmat_joined <- Qmat %>% left_join(meta, by = "Ind")
Qmat_joined <- Qmat_joined %>%
  mutate(
    POP = trimws(POP),                   # limpia espacios invisibles
    POP = iconv(POP, to = "UTF-8"),      # corrige encoding
    POP_label = recode(POP,              # estandariza nombres
                       "1" = "AA",
                       "2" = "AB",
                       "3" = "BA",
                       "4" = "BB",
                       "5" = "C",
                       "6" = "DA",
                       "7" = "DB",
                       "8" = "EA",
                       "9" = "EB",
                       "10" = "EC")
  )

# ============================================================
# 7. ASIGNAR CLUSTER MAYORITARIO
# ============================================================

Qmat_joined <- Qmat_joined %>%
  mutate(
    AssignedCluster = apply(select(., starts_with("Q")), 1, which.max)
  )

# ============================================================
# 8. CALCULAR ProbAssigned (probabilidad del cluster asignado)
# ============================================================

Qmat_joined <- Qmat_joined %>%
  mutate(
    ProbAssigned = select(., starts_with("Q"))[
      cbind(1:n(), AssignedCluster)
    ]
  )

# ============================================================
# 9. TABLA POP → CLUSTER
# ============================================================

pop_cluster <- Qmat_joined %>%
  group_by(AssignedCluster, POP_label)%>%
  summarise(
    n_ind = n(),
    meanProb = mean(ProbAssigned),
    .groups = "drop"
  ) %>%
  arrange(AssignedCluster, desc(n_ind))

cat("\n📋 Tabla POP vs Cluster:\n")
print(pop_cluster)

# ============================================================
# 10. HEATMAP POP vs CLUSTER
# ============================================================

p_heatmap <- ggplot(pop_cluster, aes(x = factor(AssignedCluster), y = factor(POP_label), fill = n_ind)) +
  geom_tile() +
  geom_text(aes(label = n_ind), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "POP vs Cluster",
       x = "Cluster asignado",
       y = "Población") +
  theme_minimal(base_size = 14)

print(p_heatmap)
dir.create(here("results","structure","plots"),showWarnings = FALSE)

ggsave(
  filename = here("results","structure","plots","heatmap_popcluster.png"),
  plot = p_heatmap,
  width = 8,
  height = 6,
  dpi = 300
)
