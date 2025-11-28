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
# CALCULAR Evanno table (versión corregida)
# ============================================================

evanno <- df %>%
  group_by(K) %>%
  summarise(
    Mean_LnP = mean(LnP),
    SD_LnP = sd(LnP)
  ) %>%
  arrange(K) %>%
  mutate(
    LnP_prime = c(NA, diff(Mean_LnP)),                    # first derivative
    LnP_doubleprime = c(NA, diff(LnP_prime[-1]), NA),     # second derivative corregido
    DeltaK = abs(LnP_doubleprime) / SD_LnP
  )

bestK <- evanno %>%
  filter(DeltaK == max(DeltaK, na.rm = TRUE)) %>%
  pull(K)

cat("\nTabla Evanno lista.\n")
print(evanno)
cat("\n")
cat("✨ Mejor K según método de Evanno:", bestK, "\n")


# ============================================================
# PALETA DE COLORES Y CONFIGURACIÓN
# ============================================================

# Paleta de colores que usamos anteriormente (azules y verdes)
mi_paleta <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Colores específicos para cada gráfico
color_lnp <- "#1f77b4"      # Azul principal
color_deltak <- "#d62728"   # Rojo para Delta K
color_boxplot <- "#2ca02c"  # Verde para boxplot
color_bestk <- "#ff7f0e"    # Naranja para línea del mejor K

# ============================================================
# 1. GRAFICO: Mean LnP(K) con barras de error (MEJORADO)
# ============================================================

p1 <- ggplot(evanno, aes(x = K, y = Mean_LnP)) +
  geom_point(size = 4, color = color_lnp, alpha = 0.8) +
  geom_line(color = color_lnp, linewidth = 1, alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_LnP - SD_LnP, ymax = Mean_LnP + SD_LnP),
                width = 0.3, color = "gray40", alpha = 0.7) +
  geom_vline(xintercept = bestK, linetype = "dashed", 
             color = color_bestk, linewidth = 1) +
  labs(
    #  title = "Mean LnP(K) con barras de error",
    #     subtitle = paste("Máxima probabilidad en K =", bestK),
    x = "Número de poblaciones (K)",
    y = "Mean Ln Probability of Data") +
  scale_x_continuous(breaks = unique(evanno$K)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.minor = element_blank()
  )

print(p1)

# ============================================================
# 2. GRAFICO: Delta K vs K (Evanno) - MEJORADO
# ============================================================

p2 <- ggplot(evanno, aes(x = K, y = DeltaK)) +
  geom_segment(aes(x = K, xend = K, y = 0, yend = DeltaK), 
               color = color_deltak, alpha = 0.6, linewidth = 1) +
  geom_point(size = 5, color = color_deltak, fill = "white", 
             shape = 21, stroke = 1.5) +
  geom_vline(xintercept = bestK, linetype = "dashed", 
             color = color_bestk, linewidth = 1.2) +
  labs(
  #    title = "Método de Evanno: ΔK",
  #     subtitle = paste("K óptimo =", bestK),
       x = "Número de poblaciones (K)",
       y = "Delta K") +
  scale_x_continuous(breaks = unique(evanno$K)) +
  annotate("label", x = bestK, y = max(evanno$DeltaK, na.rm = TRUE) * 0.9,
           label = paste("Mejor K =", bestK), 
           color = color_bestk, fill = "white", fontface = "bold") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.minor = element_blank()
  )

print(p2)

# ============================================================
# 3. GRAFICO: Distribución de LnP por K (boxplot) - MEJORADO
# ============================================================

p3 <- ggplot(df, aes(x = factor(K), y = LnP, fill = factor(K))) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, color = "gray30") +
  geom_vline(xintercept = as.numeric(factor(bestK)), 
             linetype = "dashed", color = color_bestk, linewidth = 1) +
  scale_fill_manual(values = mi_paleta, name = "K") +
  labs(
  #  title = "Distribución de LnP(K) por réplicas",
  #     subtitle = paste("Línea punteada indica K =", bestK),
       x = "Número de poblaciones (K)",
       y = "Ln Probability of Data") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p3)

# ============================================================
# 4. GRAFICO EXTRA: Comparación lado a lado (opcional)
# ============================================================

library(patchwork)

p_combined <- (p1 | p2) / p3 +
  plot_annotation(
    title = "Análisis de Estructura Poblacional - Lagartijas",
    subtitle = "Resultados del método de Evanno para determinar K óptimo",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  )

print(p_combined)

# ============================================================
# RESULTADOS NUMÉRICOS
# ============================================================

cat("\n" ,strrep("=", 60), "\n")
cat("📊 RESUMEN DE RESULTADOS - MÉTODO DE EVANNO\n")
cat(strrep("=", 60), "\n\n")

cat("K óptimo según Delta K:", bestK, "\n\n")

cat("Tabla de resultados Evanno:\n")
print(evanno)

cat("\n📈 Estadísticas por K:\n")
df %>%
  group_by(K) %>%
  summarise(
    n_replicas = n(),
    Mean_LnP = round(mean(LnP), 2),
    SD_LnP = round(sd(LnP), 2),
    DeltaK = round(first(DeltaK), 2)
  ) %>%
  print(n = Inf)

cat("\n✨ Análisis completado correctamente!\n")-

# ============================================================
# GUARDADO DE PLOTS
# ============================================================
# 
# plots_dir <- here("results", "structure", "plots")
# dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
# 
# ggsave(here("results", "structure", "plots", "plot_meanLnP.png"),
#        p1, width = 8, height = 6, dpi = 300)
# 
# ggsave(here("results", "structure", "plots", "plot_DeltaK.png"),
#        p2, width = 8, height = 6, dpi = 300)
# 
# ggsave(here("results", "structure", "plots", "plot_LnP_boxplot.png"),
#        p3, width = 8, height = 6, dpi = 300)
# 
# ggsave(here("results", "structure", "plots", "plot_combined_harvest.png"),
#        p_combined, width = 10, height = 8, dpi = 300)
# 
 cat("\n📁 Plots guardados en: results/structure/plots/\n")
