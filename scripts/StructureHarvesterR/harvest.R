library(tidyverse)
library(stringr)
library(here)

# Asume que se abre el proyecto (.Rproj) desde RStudio, por lo que getwd() ya es la raíz
results_dir <- file.path("results", "structure")

if (!dir.exists(results_dir)) {
  stop("No existe el directorio de resultados generado por run_structure_parallel.sh: ", results_dir)
}

# ============================================================
#  LEER TODOS LOS ARCHIVOS _f DE STRUCTURE
# ============================================================
files <- list.files(results_dir, pattern = "_f$", full.names = TRUE)
if (length(files) == 0) {
  stop("No se encontraron archivos _f en: ", results_dir)
}

extract_lnp <- function(path) {
  txt <- readLines(path, warn = FALSE)
  
  # linea que contiene la log-probabilidad
  lnL_line <- txt[grep("Estimated Ln Prob of Data", txt)]
  lnL_val <- as.numeric(str_extract(lnL_line, "-?\\d+\\.\\d+"))
  
  # extraer K y réplica del nombre
  fname <- basename(path)
  K <- as.numeric(str_extract(fname, "(?<=K)\\d+(?=_)"))
  rep <- as.numeric(str_extract(fname, "(?<=rep)\\d+"))
  
  tibble(
    K = K,
    Replicate = rep,
    LnP = lnL_val
  )
}

df <- map_dfr(files, extract_lnp)

print(df)
cat("\nDatos cargados correctamente.\n")
# ============================================================


# ============================================================
# CALCULAR MÉTRICAS EVANNO
# ============================================================

# Calcular manualmente para evitar problemas de longitud
evanno <- df %>%
  group_by(K) %>%
  summarise(
    Mean_LnP = mean(LnP),
    SD_LnP = sd(LnP),
    .groups = 'drop'
  ) %>%
  arrange(K)

# Calcular derivadas manualmente
n <- nrow(evanno)
evanno$LnP_prime <- c(NA, diff(evanno$Mean_LnP))
evanno$LnP_doubleprime <- c(NA, diff(evanno$LnP_prime[-1]), NA)
evanno$DeltaK <- abs(evanno$LnP_doubleprime) / evanno$SD_LnP

# Ver resultado
print(evanno)
cat("\nTabla Evanno lista.\n")
# ============================================================
bestK <- evanno %>%
  filter(DeltaK == max(DeltaK, na.rm = TRUE)) %>%
  pull(K)

# ============================================================
# 1. GRAFICO: Mean LnP(K) con barras de error
# ============================================================

p1 <- ggplot(evanno, aes(x = K, y = Mean_LnP)) +
  geom_point(size = 3, color = "blue") +
  geom_line(color = "blue") +
  geom_errorbar(aes(ymin = Mean_LnP - SD_LnP, ymax = Mean_LnP + SD_LnP),
                width = 0.2, color = "gray40") +
  labs(title = "Mean LnP(K) con barra de error",
       x = "K",
       y = "Mean Ln Probability of Data") +
  geom_vline(xintercept = bestK, linetype = "dotted", color = "red") +
  theme_minimal(base_size = 14)

print(p1)


# ============================================================
# 2. GRAFICO: Delta K vs K (criterio Evanno)
# ============================================================

p2 <- ggplot(evanno, aes(x = K, y = DeltaK)) +
  geom_point(size = 3, color = "blue") +
  geom_line(color = "blue") +
  geom_vline(xintercept = bestK, linetype = "dashed", color = "red") +
  labs(title = "Evanno ΔK",
       x = "K",
       y = "Delta K") +
  theme_minimal(base_size = 14)
  

print(p2)


# ============================================================
# 3. GRAFICO: Distribución de LnP por K (boxplot)
# ============================================================

p3 <- ggplot(df, aes(factor(K), LnP)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(alpha = 0.4, width = 0.15) +
  geom_vline(xintercept = bestK, linetype = "dashed", color = "red") +
  labs(title = "Distribución de LnP(K) por réplicas",
       x = "K",
       y = "Ln Probability of Data") +
  theme_minimal(base_size = 14)

print(p3)

cat("\n\n⭐ Gráficos generados correctamente.\n")

# ============================================================
# GUARDADO DE PLOTS
# ============================================================

plots_dir <- here("results", "structure", "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(here("results", "structure", "plots", "plot_meanLnP.png"),
       p1, width = 8, height = 6, dpi = 300)

ggsave(here("results", "structure", "plots", "plot_DeltaK.png"),
       p2, width = 8, height = 6, dpi = 300)

ggsave(here("results", "structure", "plots", "plot_LnP_boxplot.png"),
       p3, width = 8, height = 6, dpi = 300)

cat("\n📁 Plots guardados en: results/structure/plots/\n")
