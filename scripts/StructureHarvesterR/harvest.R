library(tidyverse)
library(stringr)

# ============================================================
#  LEER TODOS LOS ARCHIVOS _f DE STRUCTURE
# ============================================================
files <- list.files("results", pattern = "_f$", full.names = TRUE)

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
  geom_point(size = 3, color = "red") +
  geom_line(color = "red") +
  labs(title = "Evanno ΔK",
       x = "K",
       y = "Delta K") +
  theme_minimal(base_size = 14)



p2 <- p2 + 
  geom_vline(xintercept = bestK, linetype = "dashed", color = "black") +
  annotate("text", x = bestK, y = max(evanno$DeltaK, na.rm = TRUE),
           label = paste("Best K =", bestK), vjust = -1)

print(p2)


# ============================================================
# 3. GRAFICO: Distribución de LnP por K (boxplot)
# ============================================================

p3 <- ggplot(df, aes(factor(K), LnP)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(alpha = 0.4, width = 0.15) +
  labs(title = "Distribución de LnP(K) por réplicas",
       x = "K",
       y = "Ln Probability of Data") +
  theme_minimal(base_size = 14)

print(p3)

cat("\n\n⭐ Gráficos generados correctamente.\n")
