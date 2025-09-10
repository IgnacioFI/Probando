#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=IMT2112
# Archivo de salida
#SBATCH --output=log.out
# Cola de trabajo
#SBATCH --partition=full
# Solicitud de cpus
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1

echo "start script"
date

mpic++ -std=c++11 tarea2.cpp -o a.out
mpirun ./a.out

echo "end script"
date

