import numpy as np

#Carbono
E = 338e9 #GPa
D1 = 0.00050 #m
D2 = 0.000200 #m
gamma = 1.91 #g/cm3
gamma = gamma*1000 #kg/m3

unidad = 2.6 #metro

num_capas = 5

A = np.pi*(D1**2-D2**2)/4 #Area de la seccion tubular

#Para los paneles solares
gamma_rigido = 1.1  # Densidad del panel solar
Area_panel = unidad*unidad
masa_total = Area_panel*gamma_rigido  # Masa de un panel solar

#Para el analisi termico
alpha = 0.5e-6 #Coeficiente de dilatacion termica
deltaT = 100

#Variable para agregar barras centrales
barras_centrales = True