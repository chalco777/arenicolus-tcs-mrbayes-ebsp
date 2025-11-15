#Estimaciones de diversidad y estructura genética usando microsatélites
#Adaptación del script de Fernanda Rabanal


## Instalación y carga de paquetes

#### Para instalar los paquetes necesarios, utilizamos los siguientes comandos:

#install.packages("adegenet")
#install.packages("hierfstat")
#install.packages("pegas")
#install.packages("ade4")
#install.packages("corrplot")
#install.packages("poppr")
#install.packages("here")

#### Para cargar los paquetes instalados, escribimos:
library(adegenet)
library(here)
library(hierfstat)
library(pegas)
library(ade4)
library(corrplot)
library(poppr)
library(ggplot2)
library(tidyr)
library(dplyr)

## Preparación del conjunto de datos

### Importando nuestros datos de microsatélite a R
csvdb_path <- here("data","fst_analysis","lagartijas_usat.csv")
micros_lagartijas <- read.table(csvdb_path,sep=",",header=T)
#head(micros_lagartijas)

### ⭐⭐ NUEVO: ORDENAR POBLACIONES ALFABÉTICAMENTE ⭐⭐
# Ordenar dataframe por población
micros_lagartijas <- micros_lagartijas[order(micros_lagartijas$Population), ]
# Asegurar orden alfabético en el factor
micros_lagartijas$Population <- factor(micros_lagartijas$Population)

# Verificar el orden
cat("Poblaciones en orden alfabético:", levels(micros_lagartijas$Population), "\n")

### Convirtiendo nuestra matriz de datos a diferentes formatos de análisis...
micros_lagartijas.genind <- df2genind(X=micros_lagartijas[,c(3:ncol(micros_lagartijas))], 
                                      sep=":", ncode=NULL, 
                                      ind.names = micros_lagartijas$ID,
                                      pop = micros_lagartijas$Population, 
                                      NA.char = "NA", ploidy = 2, type = "codom")
#X=micros_lagartijas[,c(3:ncol)] es así porque se ha señalado solo las columnas donde esta la data de los locus
#ind.names es nombre de los individuos, la columna
#ploidy tambien se señala
#GENIND es información genetica por individuo

micros_lagartijas.genind

pop(micros_lagartijas.genind)

micros_lagartijas.genind@pop

#### Transformando nuestro objeto genind a un objeto hierfstat
lagartijas.hierfstat <- genind2hierfstat(micros_lagartijas.genind)
lagartijas.hierfstat

#### Transformando nuestro objeto genind a un objeto genpop
lagartijas.genpop <- genind2genpop(micros_lagartijas.genind)
lagartijas.genpop
#GENPOP es informacion genetica por población

## Estimación de los parámetros de diversidad genética

### Verificando que nuestros marcadores sean polimórficos
resumen_lagartijas <- summary(micros_lagartijas.genind)
resumen_lagartijas

resumen_lagartijas$Hobs
resumen_lagartijas$Hexp

# Extraer las heterocigosidades observada y esperada del resumen
Hobs <- resumen_lagartijas$Hobs
Hexp <- resumen_lagartijas$Hexp

# Convertir ambos a vectores numéricos con nombres de locus
Hobs <- as.numeric(Hobs)
Hexp <- as.numeric(Hexp)
Locus <- adegenet::locNames(micros_lagartijas.genind)

# Crear un data frame limpio
df_plot <- data.frame(
  Locus = Locus,
  Hobs = Hobs,
  Hexp = Hexp
) %>%
  pivot_longer(cols = c(Hobs, Hexp),
               names_to = "Tipo",
               values_to = "Valor")

# Generar el gráfico
hobs_hexp <- ggplot(df_plot, aes(x = Locus, y = Valor, fill = Tipo)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02"),
                    labels = c("Hobs", "Hexp")) +
  labs(x = "Locus", y = "Diversidad (0–1)",
       fill = "Tipo de heterocigosidad",
       title = "Comparación de Hobs y Hexp por locus") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    legend.position = "top"
  )

ggsave(here("results", "divgen_struct","hobs_hexp.png"), plot = hobs_hexp)

### Obteniendo estadísticos de resumen básicos con basic.stats de hierfstat
#### Heterocigosidad observada
estbasicos_lagartijas <- basic.stats(micros_lagartijas.genind)
estbasicos_lagartijas

estbasicos_lagartijas$Ho

colMeans(estbasicos_lagartijas$Ho, na.rm=TRUE)

barplot(colMeans(estbasicos_lagartijas$Fis, na.rm=TRUE), col ="darkseagreen2",main ="Fis")

# Extract both components
perloc <- estbasicos_lagartijas$perloc
overall <- as.data.frame(t(estbasicos_lagartijas$overall))  # convert to 1-row data frame

# Export per-locus results
write.csv(
  perloc,
  file = here("results","divgen_struct","estadisticas_perlocus.csv"),
  row.names = TRUE
)

# Export overall results
write.csv(
  overall,
  file = here("results","divgen_struct","estadisticas_overall.csv"),
  row.names = TRUE
)

### Estimando la Riqueza alélica rarefaccionada (Ar) con allelic.richness de hierfstat
Ar_lagartijas <- allelic.richness(micros_lagartijas.genind)
Ar_lagartijas

colSums(Ar_lagartijas$Ar)

# Preparar datos para ggplot
ar_data <- data.frame(
  Population = names(colSums(Ar_lagartijas$Ar)),
  Allelic_Richness = colSums(Ar_lagartijas$Ar)
)

# Ordenar por riqueza alélica
ar_data <- ar_data[order(ar_data$Allelic_Richness), ]
ar_data$Population <- factor(ar_data$Population, levels = ar_data$Population)

# Crear el plot
ar_plot <- ggplot(ar_data, aes(x = Population, y = Allelic_Richness, fill = Population)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = rep("#1b9e77", nrow(ar_data))) +
  labs(
    title = "Riqueza Alélica (Allelic Richness) por Población",
    x = "Población",
    y = "Riqueza Alélica Total",
    subtitle = "Suma de riqueza alélica rarefaccionada across 27 loci"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  # Añadir valores en las barras
  geom_text(aes(label = round(Allelic_Richness, 1)), 
            vjust = -0.5, size = 3, fontface = "bold")

# Mostrar plot
print(ar_plot)

# Guardar plot
ggsave(here("results", "divgen_struct", "allelic_richness_plot.png"), 
       plot = ar_plot, width = 10, height = 6, dpi = 300)

## Probando el equilibrio de Hardy-Weinberg
round(hw.test(micros_lagartijas.genind), digits=3)

## Evaluación de la estructura genética poblacional

### Realizando el Análisis de Componentes Principales (PCA)

#### A nivel poblacional
PCApob_lagartijas <- dudi.pca(lagartijas.genpop, scannf=FALSE, nf=2)
PCApob_lagartijas

s.label(PCApob_lagartijas$li)

#### Ahora a nivel individual
escalado_lagartijas <- scaleGen(micros_lagartijas.genind, NA.method="mean")
escalado_lagartijas

PCAind_lagartijas <- dudi.pca(escalado_lagartijas, scale=FALSE, scannf=FALSE, nf=2)
PCAind_lagartijas

##### Generando la paleta de color para nuestro gráfico
Cols<-c("yellow","chocolate","green","red4","black","purple","blue")
color_pallete_function <- colorRampPalette(colors =Cols,space = "Lab")
num_colors <- nlevels(pop(micros_lagartijas.genind))
Colores <- color_pallete_function(num_colors)

### Convertir los resultados del PCA a dataframe para ggplot
pca_df <- as.data.frame(PCAind_lagartijas$li)
pca_df$Population <- pop(micros_lagartijas.genind)  # Añadir información de poblaciones
pca_df$Individuo <- rownames(pca_df)

# Verificar nombres de las columnas (por si son diferentes)
colnames(pca_df)[1:2] <- c("Axis1", "Axis2")  # Asegurar nombres consistentes

# Calcular porcentaje de varianza explicada
var_exp <- PCAind_lagartijas$eig / sum(PCAind_lagartijas$eig) * 100

### Crear el gráfico con ggplot2 - LEYENDA FUERA DEL PLOT
pca_plot <- ggplot(pca_df, aes(x = Axis1, y = Axis2, color = Population)) +
  geom_point(size = 3, alpha = 0.8, shape = 16) +
  scale_color_manual(values = Colores, name = "Poblaciones") +
  labs(
    x = paste("Axis 1 (", round(var_exp[1], 1), "% de varianza)", sep = ""),
    y = paste("Axis 2 (", round(var_exp[2], 1), "% de varianza)", sep = ""),
    title = "Análisis de Componentes Principales (PCA)",
    subtitle = "Nivel individual - Microsatélites"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12),
    legend.position = "right",  # ✅ LEYENDA FUERA, A LA DERECHA
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", color = "gray70", linewidth = 0.3),
    plot.margin = margin(1, 1, 1, 1, "cm")  # ✅ Márgenes adecuados para leyenda
  ) +
  # Líneas de referencia en 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.6) +
  # Ajustar leyenda para hacerla más compacta si es necesario
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# Mostrar el gráfico
print(pca_plot)

# Guardar el gráfico con dimensiones adecuadas para leyenda a la derecha
ggsave(here("results", "divgen_struct", "PCA_individual_ggplot_final.png"), 
       plot = pca_plot, 
       width = 12,  # ✅ Más ancho para acomodar leyenda
       height = 8, 
       dpi = 300)

# Mensaje de confirmación
cat("✅ PCA guardado en: results/divgen_struct/PCA_individual_ggplot_final.png\n")
cat("📊 Varianza explicada - PC1:", round(var_exp[1], 1), "%, PC2:", round(var_exp[2], 1), "%\n")
cat("👥 Número de poblaciones:", nlevels(pca_df$Population), "\n")
cat("🧬 Número total de individuos:", nrow(pca_df), "\n")

### Aplicando el análisis de Fst pareadas con hierfstat
FstPareadas_lagartijas <- pairwise.neifst(lagartijas.hierfstat)
FstPareadas_lagartijas

png(here("results", "divgen_struct", "FST_pareadas_corrplot.png"), 
    width = 10, height = 8, units = "in", res = 300)
corrplot(as.matrix(FstPareadas_lagartijas), 
         is.corr = FALSE, 
         type = "lower",
         col = colorRampPalette(c("#1b9e77", "red"))(100),
         tl.col = "black", tl.srt = 45, diag = FALSE,
         cl.pos = "r", mar = c(0, 0, 1, 0))
dev.off()

### Realizando el análisis de AMOVA con poppr
Poblaciones <- data.frame("Population"=micros_lagartijas$Population)
head(Poblaciones)

strata(micros_lagartijas.genind) <- Poblaciones
setPop(micros_lagartijas.genind) <- ~Population
micros_lagartijas.genind

table(strata(micros_lagartijas.genind, ~Population))

AMOVA_lagartijas <- poppr.amova(micros_lagartijas.genind, ~Population, within = FALSE)
AMOVA_lagartijas




