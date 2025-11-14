#install.packages("here")
#install.packages("dplyr")

# ==============================
# 1. Leer archivo original (sin editar)
# ==============================
library(here)
library(dplyr)

path_db = here("data","microsatellite_genotypes.gen.txt")
lines <- readLines(path_db)

# ==============================
# 2. Inicializar variables
# ==============================
pop_num <- 0
processed_lines <- c()

# ==============================
# 3. Recorrer líneas
# ==============================
for (ln in lines) {
  # Si la línea dice "pop", incrementa población y continúa
  if (grepl("^\\s*pop\\s*$", ln, ignore.case = TRUE)) {
    pop_num <- pop_num + 1
  } else if (grepl(",", ln)) {
    # 🔹 Reemplazar SOLO la coma por el número de población (sin insertar tabs)
    ln <- sub(",", as.character(pop_num), ln)
    processed_lines <- c(processed_lines, ln)
  } else if (grepl("^ID", ln)) {
    # Mantener cabecera igual
    processed_lines <- c(processed_lines, ln)
  }
}
# ==============================
# 4. Agregar header personalizado antes de escribir
# ==============================
header_line <- paste(
  "ID", "pop", "sa5220", "sa5229", "sa6002", "sa6012", "sar80", "sar84",
  "sarms0001", "sarms0344", "sarms0473", "sarms0506", "sarms0739", "sarms0830",
  "sarms2770", "sarms3213", "sarms3490", "sarms3645", "sarms4015", "sarms4354",
  "sarms4545", "sarms4547", "sarms5185", "sarms5711", "sarms5839", "sarms5968",
  "sarms6064", "sarms6346", "sarms7111", sep = "\t"
)

processed_lines <- c(header_line, processed_lines)

# ==============================
# 5. Guardar archivo limpio
# ==============================
txtclean_path <- here("data","fst_analysis","microsatellite_genotypes.gen.txt")
writeLines(processed_lines, txtclean_path)

# ==============================
# 6. Leer el nuevo archivo como tabla
# ==============================
lagartijas_raw <- read.table(txtclean_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)

#Añadir columna de poblacion desde sample_to_region.tsv
path_regions <- here("data", "sample_to_region.tsv")
regions_df <- read.csv(path_regions, sep = "\t", header = TRUE)
head(regions_df)

df_joined <- merge(lagartijas_raw, regions_df, by.x = "ID", by.y = "sample_norm") %>%
  mutate(
    pop = as.character(pop)
  )

lagartijas_raw <- df_joined %>%
  mutate(
    pop = coalesce(region, pop)
  ) %>%
  dplyr::select(-region)

head(lagartijas_raw)

# Verificar
head(lagartijas_raw[,1:6])
unique(lagartijas_raw$pop)

# Función para convertir "109112" → "109:112" y "000000" → NA
convert_geno <- function(x) {
  x <- ifelse(x == "000000" | x == "0", NA, x)
  gsub("(...)(...)", "\\1:\\2", x)
}

# Identificar columnas de loci (todas excepto ID y pop)
loci_cols <- setdiff(names(lagartijas_raw), c("ID", "pop"))

# Aplicar conversión a todas las columnas de loci
lagartijas_fmt <- lagartijas_raw %>%
  mutate(across(all_of(loci_cols), convert_geno)) %>%
  rename(Population = pop)

csvdb_path <- here("data","fst_analysis","lagartijas_usat.csv")
# Guardar versión final lista para análisis
write.csv(lagartijas_fmt, csvdb_path, row.names = FALSE, quote = FALSE)

# Verificar
head(lagartijas_fmt[, 1:6])
