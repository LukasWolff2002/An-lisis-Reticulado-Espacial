import numpy as np
import h5py

A = np.zeros((3, 3))

for i in range(3):
    A[:, i] = i

print(A)

#hdf5 es un dato que es mas general para otras aplicaciones

#f = h5py.File('archivo.h5', 'w')

#f['A'] = A

#f.close()

#Puedo mirar los datos en el navegador
#Buscar hdf5 view

#Tambien puedo leer mediante python

f = h5py.File('archivo.h5', 'r')

A = f['A'][()] #Con este subcorchete final le pido que me entregue todos los datos

#Tambien puedo extraer solo una columna
B = f['A'][:, 1]

print(A)
print(B)

f.close()

#Es decir no cargo todos los datos, pido los datos segun variables
#Lo cual es mucho mas eficiente, ya que no cargo todos los datos en la memoria



