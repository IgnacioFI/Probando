import matplotlib.pyplot as plt
import pandas as pd

# Lee el archivo csv (ajusta el nombre si es diferente)
data = pd.read_csv('hist_convergencia_10.csv')

# Extrae las columnas de iteración y error
iterations = data['Iteracion']
errors = data['Error']

# Grafica el error vs iteraciones (en escala logarítmica)
plt.figure(figsize=(8, 5))
plt.plot(iterations, errors, marker='o', linestyle='-', label='Error por iteración')
plt.yscale('log')
plt.xlabel('Iteración')
plt.ylabel('Error (norma del residuo)')
plt.title('Convergencia del algoritmo (10 hilos)')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()
