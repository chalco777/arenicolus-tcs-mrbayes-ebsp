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


## Preparación del conjunto de datos

### Importando nuestros datos de microsatélite a R
csvdb_path <- here("data","fst_analysis","lagartijas_usat.csv")
micros_lagartijas <- read.table(csvdb_path,sep=",",header=T)
#head(micros_lagartijas)

### Convirtiendo nuestra matriz de datos a diferentes formatos de análisis y verificando la asignación de los individuos a sus respectivas poblaciones

#### Convirtiendo nuestra tabla de datos a formato genind con adegenet
micros_lagartijas.genind <- df2genind(X=micros_lagartijas[,c(3:ncol(micros_lagartijas))], sep=":", ncode=NULL, ind.names= micros_lagartijas$ID,
                            pop=micros_lagartijas$Population, NA.char="NA", ploidy=2, type="codom")
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

barplot(colSums(Ar_lagartijas$Ar), col ="powderblue",main ="Ar")

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

plot(PCAind_lagartijas$li,col=transp(Colores[pop(micros_lagartijas.genind)],0.8),pch=20,cex=2,xlab="Axis.1 XXX%",ylab="Axis.2 XXX%",bty="n",cex.axis=1.25,cex.lab=1.5)
abline(h=0,lty=3,col="gray75")
abline(v=0,lty=3,col="gray75")
legend("bottomleft",legend=levels(pop(lagartijas.genind)),xpd=TRUE,pt.cex=2,text.font=1,cex=1.25,col=Colores,pch=20)

### Aplicando el análisis de Fst pareadas con hierfstat
FstPareadas_lagartijas <- pairwise.neifst(lagartijas.hierfstat)
FstPareadas_lagartijas

corrplot(as.matrix(FstPareadas_lagartijas), is.corr=FALSE, type="lower")

### Realizando el análisis de AMOVA con poppr
Poblaciones <- data.frame("Population"=micros_lagartijas$Population)
head(Poblaciones)

strata(lagartijas.genind) <- Poblaciones
setPop(lagartijas.genind) <- ~Population
lagartijas.genind

table(strata(lagartijas.genind, ~Population))

AMOVA_lagartijas <- poppr.amova(lagartijas.genind, ~Population, within = FALSE)
AMOVA_lagartijas




