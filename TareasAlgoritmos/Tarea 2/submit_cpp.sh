#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=IMT2112
# Archivo de salida
#SBATCH --output=output_%j.txt
# Cola de trabajo
#SBATCH --partition=full
# Solicitud de cpus
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
# Informe
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ignacio.fullerton@uc.cl

echo "start script"
date

mpic++ tarea2.cpp -std=c++11
time mpirun ./a.out

echo "end script"
date

