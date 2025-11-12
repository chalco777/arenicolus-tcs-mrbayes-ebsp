

```bash
sudo apt update

# 1 Instalar OpenMPI (implementación típica de MPI)
sudo apt install openmpi-bin libopenmpi-dev

# 2 Instalar la versión MPI de MrBayes
sudo apt install mrbayes-mpi

# 3 check

which mb-mpi
mb-mpi -h   # para ver que arranca

# Run with 6 cores

mpirun -np 6 mb-mpi mtDNA_mrbayes.nex

```
