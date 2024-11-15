# Importación de librerías
import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import h5py

print('Revisando el diseño en modelo.h5 \n')

# Cargar los datos desde el archivo HDF5
f = h5py.File('CajaSatelite.h5', 'r')

nodes_tags = f['nodes_tags'][()]
nodes_xyz = f['nodes_xyz'][()]
nodes_fixities = f['nodes_fixities'][()]
elements_tags = f['elements_tags'][()]
elements_connectivities = f['elements_connectivities'][()]
elements_section_info = f['elements_section_info'][()]
panel_nodes = f['panel_nodes'][()]

f.close()

# Asegurarse de que nodes_tags sean enteros
nodes_tags = nodes_tags.astype(int)
nodes_xyz = nodes_xyz.astype(float)

# Asegurarse de que nodes_fixities sean enteros (0 o 1)
# Suponiendo que nodes_fixities es un array de la forma [[node_tag, fx, fy, fz], ...]
nodes_fixities = nodes_fixities.astype(int)

# Inicializar el modelo
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# Definir el material elástico
E_fibra_carbono = 338e9  # Módulo de Young de la fibra de carbono (Pa)
gamma_fibra_carbono = 1.91 * 1000  # kg/m³
ops.uniaxialMaterial('Elastic', 1, E_fibra_carbono)

# Crear un diccionario para almacenar las áreas y momentos de inercia de sección
section_areas = {}
section_inertias = {}

for i in range(len(elements_tags)):
    D2, espesor = elements_section_info[i]
    D1 = D2 + 2 * espesor
    area = (np.pi / 4) * (D1**2 - D2**2)
    inertia = (np.pi / 64) * (D1**4 - D2**4)
    
    # Usamos el índice como identificador de sección
    section_id = i + 1
    section_areas[section_id] = area
    section_inertias[section_id] = inertia

# Definir nodos y aplicar restricciones
for i, tag in enumerate(nodes_tags):
    x, y, z = nodes_xyz[i]
    ops.node(int(tag), x, y, z)
    
    # Aplicar restricciones si existen
    fixity_rows = nodes_fixities[nodes_fixities[:, 0] == tag]
    if fixity_rows.size > 0:
        # Si hay restricciones para este nodo
        fixity = fixity_rows[0][1:]  # [fx, fy, fz]
        fixity = [int(f) for f in fixity]  # Asegurar que son enteros
        ops.fix(int(tag), *fixity)
    else:
        # Si no hay restricciones, dejar el nodo libre
        pass

# Definir elementos (barras)
for i in range(len(elements_tags)):
    element_id = int(elements_tags[i])
    nodo_i = int(elements_connectivities[i][0])
    nodo_j = int(elements_connectivities[i][1])
    area = section_areas[i + 1]  # Asumiendo que section_id = i + 1
    ops.element('Truss', element_id, nodo_i, nodo_j, area, 1, '-rho', gamma_fibra_carbono)

# Inicializar diccionario de masas nodales
masa_nodal = {tag: 0.0 for tag in nodes_tags}

# Calcular masa de las barras y asignarla a los nodos
masa_total_barras = 0
for i in range(len(elements_tags)):
    element_id = int(elements_tags[i])
    nodo_i = int(elements_connectivities[i][0])
    nodo_j = int(elements_connectivities[i][1])
    
    idx_i = np.where(nodes_tags == nodo_i)[0][0]
    idx_j = np.where(nodes_tags == nodo_j)[0][0]
    coord_i = nodes_xyz[idx_i]
    coord_j = nodes_xyz[idx_j]
    
    longitud = np.linalg.norm(coord_j - coord_i)
    area = section_areas[i + 1]
    volumen = area * longitud
    masa_barra = volumen * gamma_fibra_carbono
    masa_total_barras += masa_barra
    
    masa_nodal[nodo_i] += masa_barra / 2
    masa_nodal[nodo_j] += masa_barra / 2

print(f'Masa total barras \t: {masa_total_barras:.2f} kg')

# Masa de los paneles
gamma_panel = 1.1  # kg/m²
area_panel_total = 0
masa_total_panel = 0

# Acumular la masa de los paneles en los nodos
for panel in panel_nodes:
    nodo1, nodo2, nodo3 = panel
    idx1 = np.where(nodes_tags == nodo1)[0][0]
    idx2 = np.where(nodes_tags == nodo2)[0][0]
    idx3 = np.where(nodes_tags == nodo3)[0][0]
    A = nodes_xyz[idx1]
    B = nodes_xyz[idx2]
    C = nodes_xyz[idx3]
    AB = B - A
    AC = C - A
    area_panel = 0.5 * np.linalg.norm(np.cross(AB, AC))
    masa_panel = area_panel * gamma_panel
    area_panel_total += area_panel
    masa_total_panel += masa_panel
    masa_nodal[nodo1] += masa_panel / 3
    masa_nodal[nodo2] += masa_panel / 3
    masa_nodal[nodo3] += masa_panel / 3

print(f'Area total de panel \t: {area_panel_total:.2f} m²')
print('Area requerida \t\t: 3333.33 m²')

if area_panel_total - 3333.333333 < 0:
    print('Suficiente area \t: False \n')
else:
    print('Suficiente area \t: True \n')

print(f'Masa total panel \t: {masa_total_panel:.2f} kg')

RME = (masa_total_barras / masa_total_panel) * 100
print(f'RME \t\t\t: {RME:.2f}%')
if RME > 20:
    print('Cumple RME \t\t: False')
else:
    print('Cumple RME \t\t: True')

# Asignar masas a los nodos
for tag in nodes_tags:
    masa = masa_nodal[tag]
    ops.mass(int(tag), masa, masa, masa)

# Funciones de análisis
def realizar_analisis_frecuencias(num_modes=1):
    """Realiza un análisis modal y retorna las frecuencias naturales."""
    eigenvalues = ops.eigen(num_modes)
    eigenfrequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
    return eigenfrequencies

def realizar_analisis_inercial(aceleraciones):
    # Realizar análisis inercial
    aceleracion_magnitud = 0.1 * 9.81
    aceleraciones = [
        [aceleracion_magnitud, 0, 0],
        [0, aceleracion_magnitud, 0],
        [0, 0, aceleracion_magnitud]
    ]

    fuerzas_maximas = {int(tag): 0.0 for tag in elements_tags}

    for case_id, acel in enumerate(aceleraciones, start=1):
        ops.timeSeries('Constant', case_id)
        ops.pattern('Plain', case_id, case_id)
        
        for tag in nodes_tags:
            masa = masa_nodal[tag]
            Fx = masa * acel[0]
            Fy = masa * acel[1]
            Fz = masa * acel[2]
            ops.load(int(tag), Fx, Fy, Fz)
        
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.algorithm('Linear')
        ops.analysis('Static')
        ops.analyze(1)
        
        for tag in elements_tags:
            element_id = int(tag)
            fuerzas = ops.eleResponse(element_id, 'axialForce')
            fuerza_axial = fuerzas[0]
            fuerzas_maximas[element_id] += fuerza_axial ** 2
        
        ops.remove('loadPattern', case_id)
        ops.wipeAnalysis()

    

    for key in fuerzas_maximas:
        fuerzas_maximas[key] = (fuerzas_maximas[key]) ** 0.5

    return fuerzas_maximas



# Realizar el análisis modal
eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)

print(f'Frecuencia Fundamental \t: {eigenfrequencies[0]:.3f} Hz \n')

# Realizar análisis inercial
aceleracion_magnitud = 0.1 * 9.81  # 0.1g
aceleraciones = [
    [aceleracion_magnitud, 0, 0],
    [0, aceleracion_magnitud, 0],
    [0, 0, aceleracion_magnitud]
]
fuerzas_maximas = realizar_analisis_inercial(aceleraciones)

print(f'La fuerza máxima experimentada por inercia es de {max(fuerzas_maximas.values()) / 1000:.2f} kN \n')





