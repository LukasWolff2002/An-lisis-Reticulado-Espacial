# Importación de librerías
import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import h5py
import sys

file = sys.argv[1]

print(f'Revisando el diseño en {file} \n')

# Cargar los datos desde el archivo HDF5
f = h5py.File(file, 'r')

nodes_tags = f['nodes_tags'][()]
nodes_xyz = f['nodes_xyz'][()]
nodes_fixities = f['nodes_fixities'][()]
elements_tags = f['elements_tags'][()]
elements_connectivities = f['elements_connectivities'][()]
elements_section_info = f['elements_section_info'][()]
panel_nodes = f['panel_nodes'][()]

f.close()

#Vamos a verificar si cumple con los requerimientos de diseño
nodos_fijos = 0
for elemenot in nodes_fixities:
    if elemenot[1] == 1 or elemenot[2] == 1 or elemenot[3] == 1:
        nodos_fijos += 1
        
nodos_apoyos = nodes_tags[8:nodos_fijos]
a = 0
for nodo in nodos_apoyos:
    
    coordenada = nodes_xyz[nodo-1]
    
    if abs(coordenada[0]) > 3.9 or abs(coordenada[1]) > 3.3 or abs(coordenada[2]) > 1.3:
        a = 1

if a != 0:
    print('Cumple soportes \t: False \n')

else:
    print('Cumple soportes \t: True \n')
    



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
masa_total_barras=0
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
print(f'Masa total estructura \t: {masa_total_barras:.2f} kg')

import numpy as np

def calcular_masa_centroide_paneles(panel_nodes, nodes_tags, nodes_xyz, gamma_panel):
    """
    Calcula la masa y el centroide (x, y, z) de cada panel.

    Parámetros:
        panel_nodes (list of list of int): Lista de paneles, cada uno definido por los tags de sus nodos [nodo1, nodo2, nodo3].
        nodes_tags (array-like of int): Lista o array de tags de nodos.
        nodes_xyz (array-like of list/tuple of float): Lista o array de coordenadas de nodos [x, y, z], alineadas con nodes_tags.
        gamma_panel (float): Densidad superficial (kg/m²) de los paneles.

    Retorna:
        list of dict: Lista donde cada diccionario contiene 'panel_id', 'masa', y 'centroide' del panel.
    """
    # Convertir nodes_tags y nodes_xyz a arrays de NumPy para facilitar la indexación
    nodes_tags = np.array(nodes_tags)
    nodes_xyz = np.array(nodes_xyz)

    panel_info = []  # Lista para almacenar la información de cada panel

    for i, panel in enumerate(panel_nodes):
        try:
            # Verificar que el panel tiene exactamente 3 nodos (triangular)
            if len(panel) != 3:
                print(f"Panel {i+1} no tiene exactamente 3 nodos. Saltando.")
                continue

            # Obtener los tags de los nodos del panel
            nodo1, nodo2, nodo3 = panel

            # Encontrar los índices de los nodos en nodes_tags
            idx1 = np.where(nodes_tags == nodo1)[0]
            idx2 = np.where(nodes_tags == nodo2)[0]
            idx3 = np.where(nodes_tags == nodo3)[0]

            # Verificar que los nodos existen
            if len(idx1) == 0 or len(idx2) == 0 or len(idx3) == 0:
                print(f"Uno o más nodos del panel {i+1} no existen. Saltando.")
                continue

            idx1 = idx1[0]
            idx2 = idx2[0]
            idx3 = idx3[0]

            # Obtener las coordenadas de los nodos
            A = nodes_xyz[idx1]
            B = nodes_xyz[idx2]
            C = nodes_xyz[idx3]

            # Calcular el centroide del panel
            centroide = (A + B + C) / 3

            # Calcular el área del triángulo utilizando el producto cruz
            AB = B - A
            AC = C - A
            area = 0.5 * np.linalg.norm(np.cross(AB, AC))

            # Calcular la masa del panel
            masa = area * gamma_panel

            # Almacenar la información en un diccionario
            panel_info.append({
                'panel_id': i+1,             # Identificador del panel (puedes modificarlo según tus necesidades)
                'masa': masa,                 # Masa del panel en kg
                'centroide': tuple(centroide)  # Coordenadas del centroide (x, y, z)
            })

        except Exception as e:
            print(f"Error al procesar el panel {i+1}: {e}")

    I_paneles = 0

    for panel in panel_info:
        I_paneles +=(panel['masa']/3)*((panel['centroide'][0]**2)+(panel['centroide'][1]**2)+(panel['centroide'][2]**2))


    return I_paneles

I_paneles = calcular_masa_centroide_paneles(panel_nodes, nodes_tags, nodes_xyz, 1.1)

def calcular_masa_centroide_barras(elements_tags, elements_connectivities, nodes_tags, nodes_xyz, section_areas, gamma_fibra_carbono):
    """
    Calcula la masa y el centroide (x, y, z) de cada barra.

    Parámetros:
        elements_tags (array-like of int): Lista de identificadores de las barras.
        elements_connectivities (array-like of list of int): Conectividades de las barras [[nodo_i, nodo_j], ...].
        nodes_tags (array-like of int): Lista de tags de nodos.
        nodes_xyz (array-like of list of float): Lista de coordenadas de nodos [x, y, z].
        section_areas (dict): Diccionario donde las claves son IDs de barras y los valores son áreas de sección (m²).
        gamma_fibra_carbono (float): Densidad de la fibra de carbono (kg/m³).

    Retorna:
        list of dict: Lista donde cada diccionario contiene 'barra_id', 'masa', y 'centroide' de cada barra.
    """
    # Convertir arrays a NumPy para facilitar la indexación
    nodes_tags = np.array(nodes_tags)
    nodes_xyz = np.array(nodes_xyz)

    barras_info = []  # Lista para almacenar la información de cada barra

    for i, barra_id in enumerate(elements_tags):
        try:
            # Obtener los nodos que conectan la barra
            nodo_i, nodo_j = elements_connectivities[i]

            # Encontrar los índices de los nodos en nodes_tags
            idx_i = np.where(nodes_tags == nodo_i)[0]
            idx_j = np.where(nodes_tags == nodo_j)[0]

            # Verificar que los nodos existen
            if len(idx_i) == 0 or len(idx_j) == 0:
                print(f"Nodo faltante en la barra {barra_id}. Saltando.")
                continue

            idx_i = idx_i[0]
            idx_j = idx_j[0]

            # Obtener las coordenadas de los nodos
            coord_i = nodes_xyz[idx_i]
            coord_j = nodes_xyz[idx_j]

            # Calcular la longitud de la barra
            longitud = np.linalg.norm(coord_j - coord_i)

            # Obtener el área de la sección transversal
            area = section_areas.get(barra_id)
            if area is None:
                raise ValueError(f"Área no encontrada para la barra {barra_id}")

            # Calcular el volumen de la barra
            volumen = area * longitud

            # Calcular la masa de la barra
            masa = volumen * gamma_fibra_carbono

            # Calcular el centroide de la barra (promedio de las coordenadas de los nodos)
            centroide = (coord_i + coord_j) / 2

            # Almacenar la información en un diccionario
            barras_info.append({
                'barra_id': barra_id,          # Identificador de la barra
                'masa': masa,                  # Masa de la barra en kg
                'centroide': tuple(centroide)  # Coordenadas del centroide (x, y, z)
            })

        except Exception as e:
            print(f"Error al procesar la barra {barra_id}: {e}")

    I_barras = 0
    for barra in barras_info:
        I_barras +=(barra['masa']/3)*((barra['centroide'][0]**2)+(barra['centroide'][1]**2)+(barra['centroide'][2]**2))
        

    return I_barras

I_barras = calcular_masa_centroide_barras(elements_tags, elements_connectivities, nodes_tags, nodes_xyz, section_areas, gamma_fibra_carbono)

print(f'Iner total panel \t: {I_paneles} kg*m2')
print(f'Iner total estructura \t: {I_barras} kg*m2')

RME = (masa_total_barras / masa_total_panel) * 100
print(f'RME \t\t\t: {RME:.2f}%')
if RME > 20:
    print('Cumple RME \t\t: False \n')
else:
    print('Cumple RME \t\t: True \n')

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

# Realizar análisis inercial
aceleracion_magnitud = 0.1 * 9.81  # 0.1g
aceleraciones = [
    [aceleracion_magnitud, 0, 0],
    [0, aceleracion_magnitud, 0],
    [0, 0, aceleracion_magnitud]
]
fuerzas_maximas = realizar_analisis_inercial(aceleraciones)

FS = 2
fluencia_fibra_carbono = 2.02e9
def revisar_falla_barras ():
    factor_utilizacion = {}
    pandeo = {}
    for i in range(len(elements_tags)):
        area = section_areas[i + 1]
        tension =  (fuerzas_maximas[i+1]*FS)/area
        FU = tension/fluencia_fibra_carbono
        factor_utilizacion[elements_tags[i]]= FU

        #Ahora reviso el pandeo
        nodo_i = int(elements_connectivities[i][0])
        nodo_j = int(elements_connectivities[i][1])

        coord_i = np.array(ops.nodeCoord(nodo_i))
        coord_j = np.array(ops.nodeCoord(nodo_j))

        vector = coord_j - coord_i
        largo_barra = np.linalg.norm(vector)

        p = (fuerzas_maximas[i+1]*FS*(largo_barra**2))/((np.pi**2)*section_inertias[i+1])

        pandeo[elements_tags[i]] = p

       

    return factor_utilizacion, pandeo





factor_utilizacion, pandeo = revisar_falla_barras()

fu_max = max(factor_utilizacion.values())
if fu_max < 1:
    print(f'Cumple resistencia \t: True 0.0 < FU < {fu_max}')

else: 
    print(f'Cumple resistencia \t: False 0.0 < FU < {fu_max}')

cumple_pandeo = True
for esfuerzos in pandeo:
    if esfuerzos > E_fibra_carbono:
        cumple_pandeo = False
    
if cumple_pandeo:
    print('Cumple pandeo \t\t: True \n')

else:
    print('Cumple pandeo \t\t: False \n')





#Analisis termico
def analisis_termico ():
    
    alpha = -0.5e-6
    delta_T = 150
    thermal_load = alpha * delta_T
    matTag_with_thermal = 3
    ops.uniaxialMaterial('InitStrainMaterial', matTag_with_thermal, 1, thermal_load)

    #Obtengo los nodos de los paneles
    nodos_paneles_set = set()
    for panel in panel_nodes:
        for nodo in panel:
            nodos_paneles_set.add(nodo)
    
    #Filtro las barras que conectan los nodos de los paneles
    barras_paneles = []
    for i in range(len(elements_tags)):
        nodo_i = int(elements_connectivities[i][0])
        nodo_j = int(elements_connectivities[i][1])
        if nodo_i in nodos_paneles_set and nodo_j in nodos_paneles_set:
            barras_paneles.append([int(elements_tags[i]), nodo_i, nodo_j])

    #Aplico la carga termica a las barras de los paneles
    barra_actual = elements_tags[-1]
    
    for i in range(len(barras_paneles)):
        
        area = section_areas[i + 1]
        barra_actual += 1
        ops.element('Truss', int(barra_actual), int(barras_paneles[i][1]), int(barras_paneles[i][2]), float(area), matTag_with_thermal, '-rho', gamma_fibra_carbono)

    #Genero el analisis
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

def calcular_vector_normal_panel(puntos):
    """Calcula el vector normal de un panel definido por tres puntos."""
    A = puntos[0]
    B = puntos[1]
    C = puntos[2]
    AB = B - A
    AC = C - A
    normal = np.cross(AB, AC)
    norma = np.linalg.norm(normal)
    if norma == 0:
        return np.array([0, 0, 0])
    normal_unitario = normal / norma
    return normal_unitario

def calcular_vectores_normales_paneles(conexiones_paneles):
    """Calcula los vectores normales antes y después de la deformación para cada panel."""
    lista_normales_antes = []
    lista_normales_despues = []

    # Obtener los tags de los nodos y sus coordenadas originales
    node_tags = ops.getNodeTags()
    coords_originales = {tag: np.array(ops.nodeCoord(tag)) for tag in node_tags}

    # Obtener los desplazamientos de los nodos
    desplazamientos = {tag: np.array(ops.nodeDisp(tag)) for tag in node_tags}

    # Calcular las coordenadas deformadas
    coords_deformadas = {tag: coords_originales[tag] + desplazamientos[tag] for tag in node_tags}

    for panel in conexiones_paneles:
        # Verificar que el panel tenga al menos 3 nodos
        if len(panel) < 3:
            print(f"El panel {panel} tiene menos de 3 nodos y no se puede calcular un vector normal.")
            continue

        # Obtener las coordenadas de los nodos antes y después de la deformación
        puntos_antes = [coords_originales[nodo] for nodo in panel]
        puntos_despues = [coords_deformadas[nodo] for nodo in panel]

        # Calcular el vector normal antes y después de la deformación
        normal_antes = calcular_vector_normal_panel(puntos_antes)
        normal_despues = calcular_vector_normal_panel(puntos_despues)

        lista_normales_antes.append(normal_antes)
        lista_normales_despues.append(normal_despues)

    return lista_normales_antes, lista_normales_despues

def calcular_angulo_entre_vectores(v1, v2):
    """Calcula el ángulo en grados entre dos vectores."""
    dot_product = np.dot(v1, v2)
    norma_v1 = np.linalg.norm(v1)
    norma_v2 = np.linalg.norm(v2)
    if norma_v1 == 0 or norma_v2 == 0:
        return 0.0
    cos_theta = dot_product / (norma_v1 * norma_v2)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    angulo_rad = np.arccos(cos_theta)
    angulo_deg = np.degrees(angulo_rad)
    return angulo_deg
        
        
    

analisis_termico()
lista_normal_inicial, lista_normal_final = calcular_vectores_normales_paneles(panel_nodes)

angulo_max = 0
for i in range(len(lista_normal_inicial)):
    angulo = calcular_angulo_entre_vectores(lista_normal_inicial[i], lista_normal_final[i])
    if angulo > angulo_max:
        angulo_max = angulo
    
print(f'Frecuencia Fundamental \t: {eigenfrequencies[0]:.3f} Hz')
if eigenfrequencies[0] > 0.1:
    print('Cumple f1 \t\t: True')
else:
    print('Cumple f1 \t\t: False')
print(f'θ_all_max \t\t: {angulo_max}')
if angulo_max <= 2:
    print('Cumple θ_max \t\t: True')

else:
    print('Cumple θ_max \t\t: False')






