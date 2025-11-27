#!/bin/bash

# Configuración para tu dataset de lagartijas
MIN_K=1
MAX_K=10
REPLICAS=5
THREADS=4  # Usar solo 4 núcleos para evitar crashes

echo "=================================================="
echo "PARALELIZACIÓN DE STRUCTURE - DATASET LAGARTIJAS"
echo "=================================================="
echo "Individuos: 237"
echo "Loci: 27"
echo "Valores de K: $MIN_K a $MAX_K"
echo "Réplicas por K: $REPLICAS"
echo "Hilos paralelos: $THREADS"
echo "Burn-in: 500,000"
echo "Repeticiones: 1,000,000"
echo "=================================================="

# Crear directorios para resultados
mkdir -p results logs

# Función para ejecutar Structure
run_structure() {
    local K=$1
    local rep=$2
    local seed=$RANDOM
    
    echo "Iniciando: K=$K, réplica $rep (semilla: $seed)"
    
    nice -n 10 ./structure \
        -m mainparams \
        -e extraparams \
        -i lagartijas_STRUCTURE_fixed.txt \
        -K "$K" \
        -D "$seed" \
        -o "results/K${K}_rep${rep}" \
        > "logs/K${K}_rep${rep}.log" 2>&1
    
    echo "✅ Completado: K=$K, réplica $rep"
}

# Exportar función para GNU parallel
export -f run_structure

# Ejecutar en paralelo
echo "Iniciando análisis paralelo (4 núcleos, seguro)..."
parallel --progress --delay 1 -j $THREADS run_structure {1} {2} ::: $(seq $MIN_K $MAX_K) ::: $(seq 1 $REPLICAS)

echo "=================================================="
echo "ANÁLISIS COMPLETADO"
echo "=================================================="
echo "Resultados en: results/"
echo "Logs en: logs/"
echo "=================================================="

# Generar resumen
echo "RESUMEN DE ARCHIVOS GENERADOS:"
ls -la results/ | grep -E "(K[0-9]_rep[0-9]|_f)$"
