library(tidyverse)

# 1. Cargar el CSV original
df <- read.csv(here("data","fst_analysis","lagartijas_usat.csv"), stringsAsFactors = FALSE)

# 2. Extraer ID y Population
ID <- df$ID
Population <- df$Population
df <- df %>% select(-ID, -Population)

# 3. Separar cada locus "A:B" → dos columnas "A" y "B"
df_separated <- df %>%
  mutate(across(everything(), ~strsplit(.x, ":", fixed = TRUE))) %>%
  unnest_wider(everything(), names_sep = "_")

# 4. MODIFICACIÓN DIRECTA DE NOMBRES
current_names <- names(df_separated)
new_names <- ifelse(str_ends(current_names, "_1"),
                    str_remove(current_names, "_1"),  # Quitar _1
                    "")                               # _2 → vacío

names(df_separated) <- new_names

# 5. Crear dataframe final
df_final <- cbind(
  ID = ID,
  Population = Population,
  df_separated
)

# 6. Exportar
output_file <- here("data","fst_analysis","lagartijas_STRUCTURE.txt")

# Escribir headers (nombres de columnas)
write.table(t(c("", "", names(df_separated))),
            file = output_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Escribir datos
write.table(df_final,
            file = output_file,
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)

cat("✅ FORMATO STRUCTURE CORREGIDO\n")
cat("📁 Archivo:", output_file, "\n")

# Verificación final
cat("🔍 VERIFICACIÓN:\n")
readLines(output_file, n = 2) %>% walk(~cat(., "\n"))