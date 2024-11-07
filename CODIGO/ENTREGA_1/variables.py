import numpy as np

#Carbono
E = 338e9 #GPa
<<<<<<< HEAD:CODIGO/variables.py
D1=5.8 #mm
D2 = 2.4   #mm
D1 = D1/1000 #m
D2 = D2/1000 #m
=======
D1 = 0.0050 #m
D2 = 0.00250 #m
>>>>>>> 1c3443cdc39cfd986cb7b6c8a99076cfc8ccb5ba:CODIGO/ENTREGA_1/variables.py
gamma = 1.91 #g/cm3
gamma = gamma*1000 #kg/m3

unidad = 2.6 #metro

num_capas = 16

A = np.pi*(D1**2-D2**2)/4 #Area de la seccion tubular

#Para los paneles solares
gamma_rigido = 1.1  # Densidad del panel solar
Area_panel = unidad*unidad
masa_total = Area_panel*gamma_rigido  # Masa de un panel solar
<<<<<<< HEAD:CODIGO/variables.py

print(A)
=======
>>>>>>> 1c3443cdc39cfd986cb7b6c8a99076cfc8ccb5ba:CODIGO/ENTREGA_1/variables.py
