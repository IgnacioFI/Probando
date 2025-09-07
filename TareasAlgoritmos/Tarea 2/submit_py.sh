#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=IMT2112
# Archivo de salida
#SBATCH --output=output_%j.txt
# Cola de trabajo
#SBATCH --partition=full
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

echo "start script"
date

which python
time python tarea2_generate_matrix.py

echo "end script"
date

