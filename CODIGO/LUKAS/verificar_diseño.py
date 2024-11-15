import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import h5py

print('Revisando el diseño en modelo.h5 \n')

f = h5py.File('CajaSatelite.h5', 'r')

nodes_tags = f['nodes_tags'][()]
nodes_xyz = f['nodes_xyz'][()]
nodes_fixities = f['nodes_fixities'][()]
elements_tags = f['elements_tags'][()]
elements_connectivities = f['elements_connectivities'][()]
elements_section_info = f['elements_section_info'][()]
panel_nodes = f['panel_nodes'][()]

#Los nodos de la caja corresponden a los 8 primeros nodos
nodos_caja = nodes_tags[:8]
coordenadas_caja = nodes_xyz[:8]
nodos_apoyos = nodes_fixities[8:]

#Debo verificar que nodos apoyos esten apoyados en alguna de las caras generadas por los nodos de la caja

coordenadas_apoyos = nodes_xyz[8:8+len(nodos_apoyos)]
verificar = 0
for i in range(len(nodos_apoyos)):
    if abs(coordenadas_apoyos[i][0]) >  3.9: #No cumple el eje x
        print('El nodo', nodos_apoyos[i], 'no cumple con la condición de apoyo en x')
        verificar += 1
    if abs(coordenadas_apoyos[i][1]) >  3.3: #No cumple el eje y
        print('El nodo', nodos_apoyos[i], 'no cumple con la condición de apoyo en y')
        verificar += 1
    if abs(coordenadas_apoyos[i][2]) >  1.3: #No cumple el eje z
        print('El nodo', nodos_apoyos[i], 'no cumple con la condición de apoyo en z')
        verificar += 1

if verificar != 0:
    print('Cumple soportes \t: False')
else:
    print('Cumple soportes \t: True')

#Ahora calculo el area total del panel
def area_tres_nodos_por_numero(nodo1, nodo2, nodo3):
    """Calcula el área definida por tres nodos."""
    A = nodes_xyz[nodo1-1]
    B = nodes_xyz[nodo2-1]
    C = nodes_xyz[nodo3-1]
    AB = B - A
    AC = C - A
    area_ABC = 0.5 * np.linalg.norm(np.cross(AB, AC))
    return area_ABC
area_panel = 0
for nodos in panel_nodes:
    area_panel += area_tres_nodos_por_numero(*nodos)

print(f'Area total de panel \t: {area_panel:.2f} m2')
print('Area requerida \t\t: 3333.333333 m2')

if area_panel - 3333.333333 < 0:
    print('Suficiente area \t: False \n')
else:
    print('Suficiente area \t: True \n')

masa_total_panel = area_panel * 1.1
print(f'Masa total panel \t: {masa_total_panel:.2f} kg')
#AHora calculo la masa total de las barras


