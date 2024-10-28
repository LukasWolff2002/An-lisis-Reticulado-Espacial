import numpy as np

#Carbono
E = 338e9 #GPa
D1=5.8 #mm
D2 = 2.4   #mm
D1 = D1/1000 #m
D2 = D2/1000 #m
gamma = 1.91 #g/cm3
gamma = gamma*1000 #kg/m3

unidad = 2.6 #metro

num_capas = 2

A = np.pi*(D1**2-D2**2)/4 #Area de la seccion tubular

#Para los paneles solares
gamma_rigido = 1.1  # Densidad del panel solar
Area_panel = unidad*unidad
masa_total = Area_panel*gamma_rigido  # Masa de un panel solar

print(A)