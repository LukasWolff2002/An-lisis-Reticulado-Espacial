import numpy as np
import openseespy.opensees as ops
import pyvista as pv
from math import atan

# Declaración de variables globales
nodo_actual, barra_actual = 0, 0
gamma_fibra_carbono = 1.91 * 1000
E_fibra_carbono = 338e9
gamma_panel = 1.1

# Dimensiones de las secciones de las barras
D1_Rope, D2_Rope, A_Rope = 0.004, 0.000, None
D1_Small, D2_Small, A_Small = 0.02, 0.0175, None
D1_Medium, D2_Medium, A_Medium = 0.04, 0.0375, None
D1_Large, D2_Large, A_Large = 0.06, 0.05, None

# Parámetros geométricos
ancho_barras = 2
alto_barras = 3.9
largo_inicial_barras = 7
largo_barras = 33.96
espaciamiento = 5.66

# Listas y diccionarios globales
conexiones_paneles = []
masa_acumulada_por_nodo = {}
barras = []

# Altura y espaciado de la torre
altura_torre = 6
espaciado_torre = 2

# Variacion Termica
alpha_T = -0.5e-6
deltaT = 150


def inicializar_modelo():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

def definir_material(material_id, E):
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)
    ops.uniaxialMaterial('MinMax', 2, 1, '-min', 0)

def definir_nodos_fijos(nodes):
    global nodo_actual
    for n in nodes:
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        ops.fix(nodo_actual, 1, 1, 1)

def definir_seccion_tubo(D1, D2):
    if D1 <= D2:
        raise ValueError("El diámetro exterior debe ser mayor que el diámetro interior.")
    A = np.pi * (D1**2 - D2**2) / 4
    return A

def angulos_barras():
    alpha_2 = 45
    alpha_4 = alpha_2 + 90
    alpha_6 = alpha_4 + 90
    alpha_8 = alpha_6 + 90
    return [alpha_2, alpha_4, alpha_6, alpha_8]

def generar_elemento_axial(nodo_1, nodo_2, Tamaño):
    global barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    area = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño]
    ops.element('Truss', barra_actual, nodo_1, nodo_2, area, 1, '-rho', gamma_fibra_carbono)
    barras.append([barra_actual, nodo_1, nodo_2, Tamaño])
    barra_actual += 1

def coordenadas_cartesianas(largo, angulo_xy, altura_z):
    angulo_rad = np.radians(angulo_xy)
    x = largo * np.cos(angulo_rad)
    y = largo * np.sin(angulo_rad)
    z = altura_z
    return x, y, z

def calcular_nuevo_punto_transversal(angulo_xy, largo, distancia_transversal):
    angulo_rad = np.radians(angulo_xy)
    x_base = largo * np.cos(angulo_rad)
    y_base = largo * np.sin(angulo_rad)
    angulo_perpendicular_rad = angulo_rad + np.pi / 2
    x_nuevo = x_base + distancia_transversal * np.cos(angulo_perpendicular_rad)
    y_nuevo = y_base + distancia_transversal * np.sin(angulo_perpendicular_rad)
    return x_nuevo, y_nuevo

def nodos_apoyos(nodos_apoyo, nodos_barra):
    global nodo_actual
    for i, n in enumerate(nodos_apoyo):
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        ops.fix(nodo_actual, 1, 1, 1)
    nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

def nodos_base_barras(nodos_iniciales_barra, nodos_barra, alpha):
    global nodo_actual, ancho_barras, alto_barras, largo_inicial_barras
    x, y, z = coordenadas_cartesianas(largo_inicial_barras, alpha, alto_barras/2)
    x_1, y_1 = calcular_nuevo_punto_transversal(alpha, largo_inicial_barras, ancho_barras)
    x_2, y_2 = calcular_nuevo_punto_transversal(alpha, largo_inicial_barras, -ancho_barras)
    nodo_actual += 1
    ops.node(nodo_actual, x, y, z)
    nodo_actual += 1
    ops.node(nodo_actual, x_1, y_1, z - alto_barras)
    nodo_actual += 1
    ops.node(nodo_actual, x_2, y_2, z - alto_barras)
    nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])
    for i in range(len(nodos_barra[-1])):
        j = (i + 1) % len(nodos_barra[-1])
        generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'S')
    if len(nodos_barra) > 1:
        conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'L')

def conectar_capas(nodos_capa_1, nodos_capa_2, Tamaño_Longitudinal, Tamaño_diagonal):
    global barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    area_longitudinal = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño_Longitudinal]
    area_diagonal = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño_diagonal]
    for i in range(len(nodos_capa_1)):
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], area_longitudinal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], Tamaño_Longitudinal])
        barra_actual += 1
    for i in range(len(nodos_capa_1)):
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], area_diagonal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], Tamaño_diagonal])
        barra_actual += 1
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], area_diagonal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], Tamaño_diagonal])
        barra_actual += 1

def alargar_barras(alpha, nodos_barra):
    global nodo_actual
    num_segmentos = int(largo_barras / espaciamiento)
    for a in range(num_segmentos):
        largo_adicional = (a + 1) * espaciamiento
        x, y, z = coordenadas_cartesianas(largo_inicial_barras + largo_adicional, alpha, 0)
        x_1, y_1 = calcular_nuevo_punto_transversal(alpha, largo_inicial_barras + largo_adicional, ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(alpha, largo_inicial_barras + largo_adicional, -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z + alto_barras / 2)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z - alto_barras / 2)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z - alto_barras / 2)
        nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])
        for i in range(len(nodos_barra[-1])):
            j = (i + 1) % len(nodos_barra[-1])
            generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'S')
        if a < (num_segmentos / 2) - 1:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'M')
        elif a < (3 * num_segmentos / 4) - 1:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'M', 'M')
        else:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'S', 'S')

def nodos_paneles(nodos_barra_A, nodos_barra_B):
    for i in range(len(nodos_barra_A) - 2):
        j = i + 1
        conexiones_paneles.append([nodos_barra_A[j][0], nodos_barra_B[j][0], nodos_barra_B[j+1][0]])
        conexiones_paneles.append([nodos_barra_B[j+1][0], nodos_barra_A[j+1][0], nodos_barra_A[j][0]])

def area_tres_nodos_por_numero(nodo1, nodo2, nodo3):
    A = np.array(ops.nodeCoord(nodo1))
    B = np.array(ops.nodeCoord(nodo2))
    C = np.array(ops.nodeCoord(nodo3))
    AB = B - A
    AC = C - A
    area_ABC = 0.5 * np.linalg.norm(np.cross(AB, AC))
    return area_ABC

def acumular_masa_barras_en_nodos():
    global barras, masa_acumulada_por_nodo, A_Rope, A_Large, A_Medium, A_Small, gamma_fibra_carbono
    todos_los_nodos = set()
    for element in barras:
        _, nodo_i, nodo_j, _ = element
        todos_los_nodos.add(nodo_i)
        todos_los_nodos.add(nodo_j)
    for nodo in todos_los_nodos:
        masa_nodo = 0
        for element in barras:
            _, nodo_i, nodo_j, Seccion = element
            if nodo == nodo_i or nodo == nodo_j:
                coord_i = np.array(ops.nodeCoord(nodo_i))
                coord_j = np.array(ops.nodeCoord(nodo_j))
                longitud = np.linalg.norm(coord_j - coord_i)
                area = {'S': A_Small, 'M': A_Medium, 'L': A_Large, 'Rope': A_Rope}[Seccion]
                masa_barra = longitud * area * gamma_fibra_carbono
                masa_nodo += masa_barra / 2
        masa_acumulada_por_nodo[nodo] = masa_nodo
    return masa_acumulada_por_nodo

def acumular_masa_paneles():
    area_total = 0
    masa_total = 0
    global masa_acumulada_por_nodo, conexiones_paneles, gamma_panel
    for nodos_conectados in conexiones_paneles:
        area_panel = area_tres_nodos_por_numero(*nodos_conectados)
        area_total += area_panel
        masa_panel = area_panel * gamma_panel
        masa_total += masa_panel
        masa_por_panel = masa_panel / 4
        for i in range(3):
            if i == 1:  # Nodo central
                masa_acumulada_por_nodo[nodos_conectados[i]] = masa_acumulada_por_nodo.get(nodos_conectados[i], 0) + (masa_por_panel * 2)
            else:
                masa_acumulada_por_nodo[nodos_conectados[i]] = masa_acumulada_por_nodo.get(nodos_conectados[i], 0) + masa_por_panel
    return area_total, masa_total

def asignar_masas_a_nodos():
    global masa_acumulada_por_nodo
    for nodo, masa in masa_acumulada_por_nodo.items():
        ops.mass(nodo, masa, masa, masa)

def realizar_analisis_frecuencias(num_modes):
    eigenvalues = ops.eigen(num_modes)
    eigenfrequencies = np.sqrt(eigenvalues) / (2 * np.pi)
    print(f'\nLa frecuencia natural más baja es: {eigenfrequencies[0]} Hz\n')
    return eigenfrequencies

def visualizar_panel_3_puntos(nodos_paneles, plotter, color='green'):
    for conexiones in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in conexiones])
        if len(puntos_panel) != 3:
            raise ValueError("Cada panel debe estar definido por exactamente 3 nodos.")
        surface = pv.PolyData(puntos_panel)
        surface.faces = np.hstack([[3], [0, 1, 2]])
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.3)

def visualizar_panel_4_puntos(nodos_paneles, plotter, color='green'):
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=1)

def nodos_paneles(nodos_barra_A, nodos_barra_B):
    for i in range(len(nodos_barra_A)-2):
        j = i + 1
        #conexiones_paneles.append([nodos_barra_A[j][0], nodos_barra_B[j][0], nodos_barra_B[j+1][0], nodos_barra_A[j+1][0]])
        conexiones_paneles.append([nodos_barra_A[j][0], nodos_barra_B[j][0], nodos_barra_B[j+1][0]])
        conexiones_paneles.append([nodos_barra_B[j+1][0], nodos_barra_A[j+1][0], nodos_barra_A[j][0]])

def visualizar_caja_satelite(caras_caja, barras, conexiones_paneles):
    node_tags = ops.getNodeTags()
    node_coords = [ops.nodeCoord(tag) for tag in node_tags]
    nodes = np.array(node_coords)
    plotter = pv.Plotter()
    points = pv.PolyData(nodes)
    labels = {i + 1: str(tag) for i, tag in enumerate(node_tags)}
    plotter.add_point_labels(points, labels, font_size=10, text_color="red", always_visible=True)
    visualizar_panel_4_puntos(caras_caja, plotter, color='blue')
    visualizar_panel_3_puntos(conexiones_paneles, plotter, color='gold')
    color_mapping = {
        'S': 'green',    # Barras pequeñas
        'M': 'orange',   # Barras medianas
        'L': 'red',      # Barras grandes
        'Rope': 'purple' # Otros tipos de barras
    }
    sizes_plotted = set()
    for barra in barras:
        _, nodo_i, nodo_j, Tamaño = barra
        nodo_inicio, nodo_fin = nodo_i - 1, nodo_j - 1
        coord_inicio = nodes[nodo_inicio]
        coord_fin = nodes[nodo_fin]
        linea = pv.Line(coord_inicio, coord_fin)
        color_barra = color_mapping.get(Tamaño, 'black')
        if Tamaño not in sizes_plotted:
            plotter.add_mesh(linea, color=color_barra, line_width=2, label=Tamaño)
            sizes_plotted.add(Tamaño)
        else:
            plotter.add_mesh(linea, color=color_barra, line_width=2)
    plotter.add_legend()
    plotter.show_axes()
    plotter.show()

def main():
    global gamma_fibra_carbono, E_fibra_carbono
    global A_Rope, A_Small, A_Medium, A_Large
    global D1_Rope, D2_Rope, D1_Small, D2_Small
    global D1_Medium, D2_Medium, D1_Large, D2_Large
    global barras, conexiones_paneles, nodo_actual, barra_actual
    global alpha_T, deltaT

    inicializar_modelo()
    definir_material(material_id=1, E=E_fibra_carbono)

    # Definir secciones
    A_Rope = definir_seccion_tubo(D1=D1_Rope, D2=D2_Rope)
    A_Small = definir_seccion_tubo(D1=D1_Small, D2=D2_Small)
    A_Medium = definir_seccion_tubo(D1=D1_Medium, D2=D2_Medium)
    A_Large = definir_seccion_tubo(D1=D1_Large, D2=D2_Large)

    # Definir nodos fijos
    nodes = np.array([
        [3.9, 3.3, 1.3],
        [-3.9, 3.3, 1.3],
        [-3.9, 3.3, -1.3],
        [3.9, 3.3, -1.3],
        [3.9, -3.3, 1.3],
        [-3.9, -3.3, 1.3],
        [-3.9, -3.3, -1.3],
        [3.9, -3.3, -1.3]
    ])
    definir_nodos_fijos(nodes)

    # Definir caras de la caja
    caras_caja = [
        [1, 4, 8, 5],
        [2, 3, 7, 6],
        [1, 2, 6, 5],
        [3, 4, 8, 7],
        [1, 2, 3, 4],
        [5, 6, 7, 8]
    ]

    # Definir nodos de apoyo
    nodes_apoyos_list = [
        np.array([[3.5, 3.3, 1], [0.3, 3.3, -1], [3.9, 0.9, -1]]),
        np.array([[-3.5, 3.3, 1], [-3.9, 0.9, -1], [-0.3, 3.3, -1]]),
        np.array([[-3.5, -3.3, 1], [-0.3, -3.3, -1], [-3.9, -0.9, -1]]),
        np.array([[3.5, -3.3, 1], [3.9, -0.9, -1], [0.3, -3.3, -1]])
    ]

    nodos_barras = [[] for _ in range(4)]
    for nodes_apoy, nodos_barra in zip(nodes_apoyos_list, nodos_barras):
        nodos_apoyos(nodes_apoy, nodos_barra)

    alphas = angulos_barras()
    for alpha, nodos_barra in zip(alphas, nodos_barras):
        nodos_base_barras([largo_inicial_barras, alpha, 0], nodos_barra, alpha)
        alargar_barras(alpha, nodos_barra)

    for i in range(4):
        nodos_paneles(nodos_barras[i], nodos_barras[(i + 1) % 4])


    conexiones_paneles.append([nodos_barras[0][1][0],nodos_barras[1][1][0],nodos_barras[2][1][0]])
    conexiones_paneles.append([nodos_barras[2][1][0],nodos_barras[3][1][0],nodos_barras[0][1][0]])

    # Acumular masas y realizar análisis modal
    acumular_masa_barras_en_nodos()
    masa_barras = sum(masa_acumulada_por_nodo.values())
    area_paneles, masa_paneles = acumular_masa_paneles()
    print(f'Los paneles tienen un área total de {area_paneles} y una masa total de {masa_paneles}')
    if area_paneles < 3333.333333:
        raise ValueError("El área total de los paneles no es suficiente.")

    print(f'La masa total de las barras es: {masa_barras}')
    print(f'La proporción es: {masa_barras / masa_paneles}')
    asignar_masas_a_nodos()

    #eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)

    #Analisis termico

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    thermal_load = alpha_T * deltaT

    # Crea un nuevo material que incluye la deformación térmica
    matTag_with_thermal = 3
    ops.uniaxialMaterial('InitStrainMaterial', matTag_with_thermal, 1, thermal_load)

    nodos_paneles_2 = set()
    for panel in conexiones_paneles:
        for nodo in panel:
            nodos_paneles_2.add(nodo)
    
    # Paso 2: Filtrar las barras que conectan nodos en 'nodos_paneles'
    barras_paneles = []
    for barra in barras:
        barra_id, nodo_i, nodo_j, tamaño = barra
        if nodo_i in nodos_paneles_2 and nodo_j in nodos_paneles_2:
            barras_paneles.append(barra)

    # Paso 3: Aplicar la carga térmica a las barras filtradas
    for barra in barras_paneles:

        if barra[-1] == 'Rope':
           
            thermal_force = thermal_load * A_Rope * E_fibra_carbono
        elif barra[-1] == 'S':
           
            thermal_force = thermal_load * A_Small * E_fibra_carbono
        elif barra[-1] == 'M':
           
            thermal_force = thermal_load * A_Medium * E_fibra_carbono
        elif barra[-1] == 'L':
            
            thermal_force = thermal_load * A_Large * E_fibra_carbono

        coord_i = np.array(ops.nodeCoord(barra[1]))
        coord_j = np.array(ops.nodeCoord(barra[2]))
        vector_unitario = (coord_j - coord_i) / np.linalg.norm(coord_j - coord_i)

        fuerza_nodo_i = -thermal_force * vector_unitario
        fuerza_nodo_j = thermal_force * vector_unitario

        if barra[-1] == 'Rope':
            area = A_Rope
        elif barra[-1] == 'S':
            area = A_Small
        elif barra[-1] == 'M':
            area = A_Medium
        elif barra[-1] == 'L':
            area = A_Large

        ops.element('Truss', barra_actual, barra[1], barra[2], area, 3, '-rho', gamma_fibra_carbono)
        barra_actual += 1


        

    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

    

    def visualizar_desplazamiento_termico(barras, conexiones_paneles, conexiones_caja, escala):
        plotter = pv.Plotter()
        nodes_tags = ops.getNodeTags()
        nodes_coords = np.array([ops.nodeCoord(tag) for tag in nodes_tags])
        deformaciones = [np.array([ops.nodeDisp(tag)[i] for tag in nodes_tags]) for i in range(3)]

        truss = pv.PolyData(nodes_coords)
        truss.lines = np.hstack([[2, e[1] - 1, e[2] - 1] for e in barras])

        truss_desplazado = truss.copy()
        for i in range(3):
            truss_desplazado.points[:, i] += deformaciones[i] * escala

        plotter.add_mesh(truss, color='blue', label="Estructura Original")
        plotter.add_mesh(truss_desplazado, color='orange', label="Estructura Desplazada Térmicamente")

        plotter.show_axes()
        plotter.show()

    def calcular_vector_normal_panel(puntos):
        
        # Usamos los tres primeros puntos para definir el plano
        A = puntos[0]
        B = puntos[1]
        C = puntos[2]

        # Vectores en el plano
        AB = B - A
        AC = C - A

        # Producto cruzado para obtener el vector normal
        normal = np.cross(AB, AC)

        # Normalizar el vector
        norma = np.linalg.norm(normal)
        if norma == 0:
            return np.array([0, 0, 0])
        normal_unitario = normal / norma

        return normal_unitario

    def calcular_vectores_normales_paneles(conexiones_paneles):
  
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

            # Obtener las coordenadas de los nodos antes de la deformación
            puntos_antes = [coords_originales[nodo] for nodo in panel]
            # Obtener las coordenadas de los nodos después de la deformación
            puntos_despues = [coords_deformadas[nodo] for nodo in panel]

            # Calcular el vector normal antes de la deformación
            normal_antes = calcular_vector_normal_panel(puntos_antes)
            # Calcular el vector normal después de la deformación
            normal_despues = calcular_vector_normal_panel(puntos_despues)

            lista_normales_antes.append(normal_antes)
            lista_normales_despues.append(normal_despues)

        return lista_normales_antes, lista_normales_despues

    def calcular_angulo_entre_vectores(v1, v2):
        """
        Calcula el ángulo en grados entre dos vectores.

        Parámetros:
        - v1, v2: Vectores a comparar.

        Retorna:
        - Ángulo en grados entre v1 y v2.
        """
        # Producto punto y normas
        dot_product = np.dot(v1, v2)
        norma_v1 = np.linalg.norm(v1)
        norma_v2 = np.linalg.norm(v2)
        if norma_v1 == 0 or norma_v2 == 0:
            return 0.0
        # Cálculo del ángulo
        cos_theta = dot_product / (norma_v1 * norma_v2)
        # Asegurar que el valor esté en el rango [-1, 1]
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        angulo_rad = np.arccos(cos_theta)
        angulo_deg = np.degrees(angulo_rad)
        return angulo_deg

    



        

    visualizar_desplazamiento_termico(barras, conexiones_paneles, caras_caja, 1)

    #visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

    lista_normal_inicial, lista_normal_final = calcular_vectores_normales_paneles(conexiones_paneles)
    
    for i in range(len(lista_normal_inicial)):
        angulo = calcular_angulo_entre_vectores(lista_normal_inicial[i], lista_normal_final[i])
        if angulo > 2:
            raise ValueError(f"El ángulo entre los vectores normales del panel es mayor a 2 grados")

   

    


if __name__ == "__main__":
    main()
