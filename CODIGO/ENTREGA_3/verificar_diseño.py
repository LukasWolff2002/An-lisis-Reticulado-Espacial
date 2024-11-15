# -----------------------------------
# 1. Importación de Librerías
# -----------------------------------
import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import h5py
import sys
import os

# -----------------------------------
# 2. Definición de Funciones
# -----------------------------------

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
                'panel_id': i+1,               # Identificador del panel
                'masa': masa,                   # Masa del panel en kg
                'centroide': tuple(centroide)   # Coordenadas del centroide (x, y, z)
            })

        except Exception as e:
            print(f"Error al procesar el panel {i+1}: {e}")

    I_paneles = 0

    for panel in panel_info:
        I_paneles += (panel['masa'] / 3) * (
            (panel['centroide'][0] ** 2) +
            (panel['centroide'][1] ** 2) +
            (panel['centroide'][2] ** 2)
        )

    return I_paneles

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
            # Verificar que la barra tiene exactamente 2 nodos (Truss)
            if len(elements_connectivities[i]) != 2:
                print(f"Barra {barra_id} no tiene exactamente 2 nodos. Saltando.")
                continue

            # Obtener los tags de los nodos de la barra
            nodo_i, nodo_j = elements_connectivities[i]

            # Encontrar los índices de los nodos en nodes_tags
            idx_i = np.where(nodes_tags == nodo_i)[0]
            idx_j = np.where(nodes_tags == nodo_j)[0]

            # Verificar que los nodos existen
            if len(idx_i) == 0 or len(idx_j) == 0:
                print(f"Uno o ambos nodos de la barra {barra_id} no existen. Saltando.")
                continue

            idx_i = idx_i[0]
            idx_j = idx_j[0]

            # Obtener las coordenadas de los nodos
            coord_i = np.array(nodes_xyz[idx_i])
            coord_j = np.array(nodes_xyz[idx_j])

            # Calcular la longitud de la barra
            vector = coord_j - coord_i
            largo_barra = np.linalg.norm(vector)

            # Obtener el área de la sección transversal
            area = section_areas.get(barra_id)
            if area is None:
                raise ValueError(f"Área no encontrada para la barra {barra_id}")

            # Calcular el volumen de la barra
            volumen = area * largo_barra

            # Calcular la masa de la barra
            masa = volumen * gamma_fibra_carbono

            # Calcular el centroide de la barra (punto medio)
            centroide = (coord_i + coord_j) / 2

            # Almacenar la información en un diccionario
            barras_info.append({
                'barra_id': barra_id,            # Identificador de la barra
                'masa': masa,                     # Masa de la barra en kg
                'centroide': tuple(centroide)     # Coordenadas del centroide (x, y, z)
            })

        except Exception as e:
            print(f"Error al procesar la barra {barra_id}: {e}")

    I_barras = 0
    for barra in barras_info:
        I_barras += (barra['masa'] / 3) * (
            (barra['centroide'][0] ** 2) +
            (barra['centroide'][1] ** 2) +
            (barra['centroide'][2] ** 2)
        )

    return I_barras

def realizar_analisis_frecuencias(num_modes=1):
    """Realiza un análisis modal y retorna las frecuencias naturales."""
    try:
        eigenvalues = ops.eigen(num_modes)
        eigenfrequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
        return eigenfrequencies
    except Exception as e:
        print("Error durante el análisis modal:", e)
        return []

def realizar_analisis_inercial(aceleraciones):
    """
    Realiza el análisis inercial para las aceleraciones dadas.

    Parámetros:
        aceleraciones (list of list of float): Lista de aceleraciones [ax, ay, az].

    Retorna:
        dict: Fuerzas máximas experimentadas por inercia para cada barra.
    """
    fuerzas_maximas = {int(tag): 0.0 for tag in elements_tags}

    for case_id, acel in enumerate(aceleraciones, start=1):
        try:
            # Definir la serie de tiempo y el patrón de carga
            ops.timeSeries('Constant', case_id)
            ops.pattern('Plain', case_id, case_id)

            # Aplicar las cargas inerciales a los nodos
            for tag in nodes_tags:
                masa = masa_nodal[tag]
                Fx = masa * acel[0]
                Fy = masa * acel[1]
                Fz = masa * acel[2]
                ops.load(int(tag), Fx, Fy, Fz)

            # Configurar y ejecutar el análisis
            ops.system('BandSPD')
            ops.numberer('RCM')
            ops.constraints('Plain')
            ops.integrator('LoadControl', 1.0)
            ops.algorithm('Linear')
            ops.analysis('Static')
            success = ops.analyze(1)
            if not success:
                raise RuntimeError(f'Análisis inercial fallido en el caso {case_id}')

            # Obtener las fuerzas internas en las barras
            for tag in elements_tags:
                element_id = int(tag)
                try:
                    fuerzas = ops.eleResponse(element_id, 'axialForce')
                    fuerza_axial = fuerzas[0]
                    fuerzas_maximas[element_id] += fuerza_axial ** 2
                except Exception as e:
                    print(f'Error al obtener la respuesta del elemento {element_id}:', e)

            # Limpiar cargas y análisis para la siguiente iteración
            ops.remove('loadPattern', case_id)
            ops.wipeAnalysis()

        except Exception as e:
            print(f'Error en el caso de aceleración {case_id}:', e)
            raise RuntimeError(f'Análisis inercial fallido en el caso {case_id}') from e

    # Calcular la raíz cuadrada de la suma de los cuadrados (SRSS)
    for key in fuerzas_maximas:
        fuerzas_maximas[key] = (fuerzas_maximas[key]) ** 0.5

    return fuerzas_maximas

def revisar_falla_barras():
    """
    Revisa las fallas por resistencia y pandeo en las barras.

    Retorna:
        tuple: Factor de utilización y factores de pandeo para cada barra.
    """
    factor_utilizacion = {}
    pandeo = {}
    for i in range(len(elements_tags)):
        barra_id = int(elements_tags[i])
        area = section_areas.get(barra_id, None)
        inertia = section_inertias.get(barra_id, None)
        if area is None or inertia is None:
            print(f"Área o inercia no encontrada para la barra {barra_id}. Saltando.")
            continue

        tension = (fuerzas_maximas[barra_id] * FS) / area
        FU = tension / fluencia_fibra_carbono
        factor_utilizacion[barra_id] = FU

        # Revisar pandeo
        nodo_i = int(elements_connectivities[i][0])
        nodo_j = int(elements_connectivities[i][1])

        coord_i = np.array(ops.nodeCoord(nodo_i))
        coord_j = np.array(ops.nodeCoord(nodo_j))

        vector = coord_j - coord_i
        largo_barra = np.linalg.norm(vector)

        p = (fuerzas_maximas[barra_id] * FS * (largo_barra ** 2)) / ((np.pi ** 2) * inertia)
        pandeo[barra_id] = p

    return factor_utilizacion, pandeo

def analisis_termico():
    """Realiza el análisis térmico aplicando cargas térmicas como fuerzas equivalentes nodales."""
    alpha = -0.5e-6          # Coeficiente de expansión térmica (1/K)
    delta_T = 150             # Cambio de temperatura (°C)

    # Crear un nuevo load pattern para las cargas térmicas
    thermal_load_id = len(aceleraciones) + 1  # Asegúrate de que este ID no colisione con otros
    ops.timeSeries('Constant', thermal_load_id)
    ops.pattern('Plain', thermal_load_id, thermal_load_id)

    # Obtener los nodos de los paneles
    nodos_paneles_set = set()
    for panel in panel_nodes:
        for nodo in panel:
            nodos_paneles_set.add(nodo)

    # Filtrar las barras que conectan nodos en 'nodos_paneles'
    barras_paneles = []
    for i in range(len(elements_tags)):
        nodo_i = int(elements_connectivities[i][0])
        nodo_j = int(elements_connectivities[i][1])
        if nodo_i in nodos_paneles_set and nodo_j in nodos_paneles_set:
            barras_paneles.append([int(elements_tags[i]), nodo_i, nodo_j])

    # Aplicar las cargas térmicas como fuerzas equivalentes nodales
    for barra in barras_paneles:
        element_id, nodo_i, nodo_j = barra
        area = section_areas.get(element_id, None)
        inertia = section_inertias.get(element_id, None)
        if area is None or inertia is None:
            print(f"No se encontró el área o el momento de inercia para el elemento {element_id}.")
            continue

        # Calcular la fuerza térmica
        F_T = E_fibra_carbono * area * alpha * delta_T  # Fuerza en N

        # Obtener las coordenadas de los nodos
        try:
            coord_i = np.array(ops.nodeCoord(nodo_i))
            coord_j = np.array(ops.nodeCoord(nodo_j))
        except Exception as e:
            print(f"Error al obtener coordenadas de nodos {nodo_i} o {nodo_j}: {e}")
            continue

        # Calcular la dirección unitaria del elemento
        direccion = coord_j - coord_i
        longitud = np.linalg.norm(direccion)
        if longitud == 0:
            print(f'Elemento {element_id} tiene longitud cero. Saltando.')
            continue
        direccion_unitaria = direccion / longitud

        # Descomponer la fuerza en componentes
        Fx = F_T * direccion_unitaria[0]
        Fy = F_T * direccion_unitaria[1]
        Fz = F_T * direccion_unitaria[2]

        # Aplicar las fuerzas equivalentes en los nodos
        # Nodo i recibe F_T en la dirección del elemento
        # Nodo j recibe -F_T en la dirección del elemento
        ops.load(int(nodo_i), Fx, Fy, Fz)
        ops.load(int(nodo_j), -Fx, -Fy, -Fz)

    # Configurar y ejecutar el análisis
    try:
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.algorithm('Linear')
        ops.analysis('Static')
        success = ops.analyze(1)
        if not success:
            raise RuntimeError('Análisis térmico fallido')
        print('Análisis térmico realizado exitosamente.')
    except Exception as e:
        print("Error durante el análisis térmico:", e)
        raise RuntimeError('Análisis térmico fallido') from e

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

# -----------------------------------
# 3. Ejecución Principal
# -----------------------------------
if __name__ == "__main__":
    # -----------------------------------
    # 3.1. Cargar Datos desde Archivo HDF5
    # -----------------------------------
    file = sys.argv[1]
    print(f'Revisando el diseño en {file} \n')

    # Verificar si la carpeta de destino existe, si no, crearla
    carpeta_destino = 'ENTREGA_3'
    if not os.path.exists(carpeta_destino):
        os.makedirs(carpeta_destino)

    # Cargar los datos desde el archivo HDF5
    try:
        with h5py.File(file, 'r') as f:
            nodes_tags = f['nodes_tags'][()]
            nodes_xyz = f['nodes_xyz'][()]
            nodes_fixities = f['nodes_fixities'][()]
            elements_tags = f['elements_tags'][()]
            elements_connectivities = f['elements_connectivities'][()]
            elements_section_info = f['elements_section_info'][()]
            panel_nodes = f['panel_nodes'][()]
    except Exception as e:
        print(f"Error al cargar el archivo HDF5: {e}")
        sys.exit(1)

    # -----------------------------------
    # 3.2. Verificar Requerimientos de Diseño
    # -----------------------------------
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

    # -----------------------------------
    # 3.3. Preparar Datos y Inicializar Modelo
    # -----------------------------------
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
    E_fibra_carbono = 338e9          # Módulo de Young de la fibra de carbono (Pa)
    gamma_fibra_carbono = 1.91e3     # Densidad de la fibra de carbono (kg/m³)
    fluencia_fibra_carbono = 2.02e9  # Fluencia de la fibra de carbono (Pa)
    FS = 2                            # Factor de seguridad

    ops.uniaxialMaterial('Elastic', 1, E_fibra_carbono)

    # Crear diccionarios para almacenar áreas y momentos de inercia de secciones
    section_areas = {}
    section_inertias = {}

    for i in range(len(elements_tags)):
        barra_id = int(elements_tags[i])
        D2, espesor = elements_section_info[i]
        D1 = D2 + 2 * espesor
        area = (np.pi / 4) * (D1**2 - D2**2)
        inertia = (np.pi / 64) * (D1**4 - D2**4)

        section_areas[barra_id] = area
        section_inertias[barra_id] = inertia

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
        barra_id = int(elements_tags[i])
        nodo_i = int(elements_connectivities[i][0])
        nodo_j = int(elements_connectivities[i][1])
        area = section_areas.get(barra_id, None)
        if area is None:
            print(f"Área no encontrada para la barra {barra_id}. Saltando.")
            continue
        ops.element('Truss', barra_id, nodo_i, nodo_j, area, 1, '-rho', gamma_fibra_carbono)

    # -----------------------------------
    # 3.4. Calcular Masas Nodales
    # -----------------------------------
    # Inicializar diccionario de masas nodales
    masa_nodal = {tag: 0.0 for tag in nodes_tags}

    # Calcular masa de las barras y asignarla a los nodos
    masa_total_barras = 0
    for i in range(len(elements_tags)):
        barra_id = int(elements_tags[i])
        nodo_i = int(elements_connectivities[i][0])
        nodo_j = int(elements_connectivities[i][1])

        idx_i = np.where(nodes_tags == nodo_i)[0]
        idx_j = np.where(nodes_tags == nodo_j)[0]

        if len(idx_i) == 0 or len(idx_j) == 0:
            print(f"Nodos {nodo_i} o {nodo_j} no existen. Saltando barra {barra_id}.")
            continue

        idx_i = idx_i[0]
        idx_j = idx_j[0]
        coord_i = nodes_xyz[idx_i]
        coord_j = nodes_xyz[idx_j]

        longitud = np.linalg.norm(coord_j - coord_i)
        area = section_areas.get(barra_id, None)
        if area is None:
            print(f"Área no encontrada para la barra {barra_id}. Saltando.")
            continue
        volumen = area * longitud
        masa_barra = volumen * gamma_fibra_carbono
        masa_total_barras += masa_barra

        masa_nodal[nodo_i] += masa_barra / 2
        masa_nodal[nodo_j] += masa_barra / 2

    # -----------------------------------
    # 3.5. Calcular Masas de Paneles
    # -----------------------------------
    gamma_panel = 1.1  # kg/m²
    area_panel_total = 0
    masa_total_panel = 0

    # Acumular la masa de los paneles en los nodos
    for panel in panel_nodes:
        nodo1, nodo2, nodo3 = panel
        try:
            idx1 = np.where(nodes_tags == nodo1)[0][0]
            idx2 = np.where(nodes_tags == nodo2)[0][0]
            idx3 = np.where(nodes_tags == nodo3)[0][0]
        except IndexError:
            print(f"Uno de los nodos del panel {panel} no existe. Saltando.")
            continue

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

    if (area_panel_total - 3333.333333) < 0:
        print('Suficiente area \t: False \n')
    else:
        print('Suficiente area \t: True \n')

    print(f'Masa total panel \t: {masa_total_panel:.2f} kg')
    print(f'Masa total estructura \t: {masa_total_barras:.2f} kg')

    # -----------------------------------
    # 3.6. Calcular Momento de Inercia
    # -----------------------------------
    I_paneles = calcular_masa_centroide_paneles(panel_nodes, nodes_tags, nodes_xyz, gamma_panel)

    I_barras = calcular_masa_centroide_barras(
        elements_tags,
        elements_connectivities,
        nodes_tags,
        nodes_xyz,
        section_areas,
        gamma_fibra_carbono
    )

    print(f'Iner total panel \t: {I_paneles:.2f} kg*m²')
    print(f'Iner total estructura \t: {I_barras:.2f} kg*m²')

    # -----------------------------------
    # 3.7. Calcular y Revisar RME
    # -----------------------------------
    RME = (masa_total_barras / masa_total_panel) * 100
    print(f'RME \t\t\t: {RME:.2f}%')
    if RME > 20:
        print('Cumple RME \t\t: False\n')
    else:
        print('Cumple RME \t\t: True\n')

    # -----------------------------------
    # 3.8. Asignar Masas a los Nodos
    # -----------------------------------
    for tag in nodes_tags:
        masa = masa_nodal[tag]
        ops.mass(int(tag), masa, masa, masa)  # Asignar masa en x, y, z

    # -----------------------------------
    # 3.9. Definir Funciones de Análisis
    # -----------------------------------
    # (Las funciones ya están definidas anteriormente)

    # -----------------------------------
    # 3.10. Realizar Análisis Modal
    # -----------------------------------
    eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)

    # -----------------------------------
    # 3.11. Realizar Análisis Inercial
    # -----------------------------------
    aceleracion_magnitud = 0.1 * 9.81  # 0.1g
    aceleraciones = [
        [aceleracion_magnitud, 0, 0],
        [0, aceleracion_magnitud, 0],
        [0, 0, aceleracion_magnitud]
    ]
    fuerzas_maximas = realizar_analisis_inercial(aceleraciones)

    # -----------------------------------
    # 3.12. Revisar Falla en Barras
    # -----------------------------------
    factor_utilizacion, pandeo = revisar_falla_barras()

    fu_max = max(factor_utilizacion.values()) if factor_utilizacion else 0
    if fu_max < 1:
        print(f'Cumple resistencia \t: True 0.0 < FU < {fu_max:.2f}')
    else:
        print(f'Cumple resistencia \t: False 0.0 < FU < {fu_max:.2f}')

    cumple_pandeo = True
    for barra_id, p in pandeo.items():
        if p > E_fibra_carbono:
            cumple_pandeo = False
            break

    if cumple_pandeo:
        print('Cumple pandeo \t\t: True \n')
    else:
        print('Cumple pandeo \t\t: False \n')

    # -----------------------------------
    # 3.13. Realizar Análisis Térmico
    # -----------------------------------
    analisis_termico()

    # -----------------------------------
    # 3.14. Calcular y Revisar Angulo Máximo
    # -----------------------------------
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
    print(f'θ_all_max \t\t: {angulo_max:.2f}°')
    if angulo_max <= 2:
        print('Cumple θ_max \t\t: True')
    else:
        print('Cumple θ_max \t\t: False')

    # -----------------------------------
    # 3.15. Guardar Resultados en HDF5
    # -----------------------------------
    try:
        ruta_completa = os.path.join(carpeta_destino, 'CajaSatelite_Final.h5')
        with h5py.File(ruta_completa, 'w') as hf5:
            hf5.create_dataset('nodes_tags', data=nodes_tags)
            hf5.create_dataset('nodes_xyz', data=nodes_xyz)
            hf5.create_dataset('nodes_fixities', data=nodes_fixities)
            hf5.create_dataset('elements_tags', data=elements_tags)
            hf5.create_dataset('elements_connectivities', data=elements_connectivities)
            hf5.create_dataset('elements_section_info', data=elements_section_info)
            hf5.create_dataset('panel_nodes', data=panel_nodes)
            hf5.create_dataset('frequencies', data=eigenfrequencies)
            hf5.create_dataset('fuerzas_maximas', data=list(fuerzas_maximas.values()))
            hf5.create_dataset('section_areas', data=list(section_areas.values()))
            hf5.create_dataset('section_inertias', data=list(section_inertias.values()))
            # Opcional: guardar la información de masa y centroide como strings
            hf5.create_dataset('info_paneles', data=np.string_(str(I_paneles)))
            hf5.create_dataset('info_barras', data=np.string_(str(I_barras)))
        print(f"Datos exportados exitosamente a '{ruta_completa}'\n")
    except Exception as e:
        print(f"Error al exportar los datos a '{ruta_completa}': {e}\n")
