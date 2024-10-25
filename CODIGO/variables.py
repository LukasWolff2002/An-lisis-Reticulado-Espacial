import numpy as np

#Carbono
E = 338e9 #
D1 = 150 #mm
D2 = 100 #mm
gamma = 1.91 #g/cm3
gamma = gamma*1000 #kg/m3

unidad = 2.6 #metro

num_capas = 10

A = np.pi*(D1**2-D2**2)/4 #Area de la seccion tubular

#Para los paneles solares
gamma_rigido = 1.1  # Densidad del panel solar
Area = unidad*unidad
masa_total = Area/gamma  # Masa de un panel solar

print(A)