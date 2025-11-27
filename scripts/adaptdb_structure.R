library(tidyverse)
library(here)

# 1. Cargar el CSV original
df <- read.csv(here("data","fst_analysis","lagartijas_usat.csv"), stringsAsFactors = FALSE)

# 2. Extraer ID y Population
ID <- df$ID
Population <- df$Population
df <- df %>% select(-ID, -Population)

# 3. Separar cada locus "A:B" → dos columnas
df_separated <- df %>%
  mutate(across(everything(), ~strsplit(.x, ":", fixed = TRUE))) %>%
  unnest_wider(everything(), names_sep = "_")

# 4. Renombrar columnas para formato STRUCTURE
current_names <- names(df_separated)
new_names <- ifelse(str_ends(current_names, "_1"),
                    str_remove(current_names, "_1"),  # columna _1 → nombre del locus
                    "")                               # columna _2 → nombre vacío

names(df_separated) <- new_names

# 5. Crear dataframe final
df_final <- cbind(
  ID = ID,
  Population = Population,
  df_separated
)

# 6. Exportar SIN HEADERS (porque STRUCTURE no los usa)
output_file <- here("data","fst_analysis","lagartijas_STRUCTURE.txt")

write.table(df_final,
            file = output_file,
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("✅ Archivo STRUCTURE generado SIN primera línea de headers.\n")
cat("📁 Archivo:", output_file, "\n")

# Verificación
cat("🔍 Primeras líneas:\n")
readLines(output_file, n = 2) %>% walk(~cat(., "\n"))
