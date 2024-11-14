# Importación de librerías
import numpy as np
import openseespy.opensees as ops
import pyvista as pv

# Declaración de variables globales
nodo_actual, barra_actual = 0, 0
gamma_fibra_carbono = 1.91 * 1000     # Densidad de la fibra de carbono (kg/m³)
E_fibra_carbono = 338e9               # Módulo de Young de la fibra de carbono (Pa)
gamma_panel = 1.1                     # Densidad del panel (kg/m²)

# Dimensiones de las secciones de las barras
D1_Rope, D2_Rope = 0.004, 0.000
D1_Small, D2_Small = 0.02, 0.0175
D1_Medium, D2_Medium = 0.04, 0.0375
D1_Large, D2_Large = 0.06, 0.05
A_Rope, A_Small, A_Medium, A_Large = None, None, None, None

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

# Variación Térmica
alpha_T = -0.5e-6  # Coeficiente de expansión térmica (1/°C)
deltaT = 150       # Cambio de temperatura (°C)

# Funciones de inicialización y definición de materiales
def inicializar_modelo():
    """Inicializa el modelo de OpenSees."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

def definir_material(material_id, E):
    """Define un material elástico en OpenSees."""
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)

def definir_seccion_tubo(D1, D2):
    """Calcula el área de una sección tubular."""
    if D1 <= D2:
        raise ValueError("El diámetro exterior debe ser mayor que el diámetro interior.")
    A = np.pi * (D1**2 - D2**2) / 4
    return A

# Funciones de generación de nodos y elementos
def definir_nodos_fijos(nodes):
    """Define nodos fijos en el modelo."""
    global nodo_actual
    for n in nodes:
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        ops.fix(nodo_actual, 1, 1, 1)

def generar_elemento_axial(nodo_1, nodo_2, Tamaño):
    """Genera un elemento axial (barra) entre dos nodos."""
    global barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    area = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño]
    ops.element('Truss', barra_actual, nodo_1, nodo_2, area, 1, '-rho', gamma_fibra_carbono)
    barras.append([barra_actual, nodo_1, nodo_2, Tamaño])
    barra_actual += 1

# Funciones geométricas y de cálculo de coordenadas
def coordenadas_cartesianas(largo, angulo_xy, altura_z):
    """Calcula las coordenadas cartesianas dado un largo, ángulo y altura."""
    angulo_rad = np.radians(angulo_xy)
    x = largo * np.cos(angulo_rad)
    y = largo * np.sin(angulo_rad)
    z = altura_z
    return x, y, z

def calcular_nuevo_punto_transversal(angulo_xy, largo, distancia_transversal):
    """Calcula un punto desplazado transversalmente desde una posición base."""
    angulo_rad = np.radians(angulo_xy)
    x_base = largo * np.cos(angulo_rad)
    y_base = largo * np.sin(angulo_rad)
    angulo_perpendicular_rad = angulo_rad + np.pi / 2
    x_nuevo = x_base + distancia_transversal * np.cos(angulo_perpendicular_rad)
    y_nuevo = y_base + distancia_transversal * np.sin(angulo_perpendicular_rad)
    return x_nuevo, y_nuevo

# Funciones para definir apoyos y nodos base
def nodos_apoyos(nodos_apoyo, nodos_barra):
    """Define los nodos de apoyo y los agrega a la lista de nodos de la barra."""
    global nodo_actual
    for n in nodos_apoyo:
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        ops.fix(nodo_actual, 1, 1, 1)
    nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

def nodos_base_barras(nodos_iniciales_barra, nodos_barra, alpha):
    """Define los nodos base de las barras y genera los elementos iniciales."""
    global nodo_actual, ancho_barras, alto_barras, largo_inicial_barras
    x, y, z = coordenadas_cartesianas(largo_inicial_barras, alpha, alto_barras / 2)
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

# Funciones para conectar capas y alargar barras
def conectar_capas(nodos_capa_1, nodos_capa_2, Tamaño_Longitudinal, Tamaño_diagonal):
    """Conecta dos capas de nodos con elementos longitudinales y diagonales."""
    global barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    area_longitudinal = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño_Longitudinal]
    area_diagonal = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño_diagonal]
    for i in range(len(nodos_capa_1)):
        # Elementos longitudinales
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], area_longitudinal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], Tamaño_Longitudinal])
        barra_actual += 1
        # Elementos diagonales
        indices = [((i + 1) % len(nodos_capa_2)), ((i - 1) % len(nodos_capa_2))]
        for idx in indices:
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[idx], area_diagonal, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[idx], Tamaño_diagonal])
            barra_actual += 1

def alargar_barras(alpha, nodos_barra):
    """Extiende las barras a lo largo de su longitud y genera los elementos correspondientes."""
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
        # Generar elementos axiales
        for i in range(len(nodos_barra[-1])):
            j = (i + 1) % len(nodos_barra[-1])
            generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'S')
        # Conectar capas según la posición
        if a < (num_segmentos / 2) - 1:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'M')
        elif a < (3 * num_segmentos / 4) - 1:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'M', 'M')
        else:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'S', 'S')

# Funciones para definir paneles y calcular áreas
def nodos_paneles(nodos_barra_A, nodos_barra_B):
    """Define los paneles conectando nodos entre dos barras."""
    for i in range(len(nodos_barra_A) - 2):
        j = i + 1
        conexiones_paneles.append([nodos_barra_A[j][0], nodos_barra_B[j][0], nodos_barra_B[j + 1][0]])
        conexiones_paneles.append([nodos_barra_B[j + 1][0], nodos_barra_A[j + 1][0], nodos_barra_A[j][0]])

def area_tres_nodos_por_numero(nodo1, nodo2, nodo3):
    """Calcula el área definida por tres nodos."""
    A = np.array(ops.nodeCoord(nodo1))
    B = np.array(ops.nodeCoord(nodo2))
    C = np.array(ops.nodeCoord(nodo3))
    AB = B - A
    AC = C - A
    area_ABC = 0.5 * np.linalg.norm(np.cross(AB, AC))
    return area_ABC

# Funciones para acumulación de masas
def acumular_masa_barras_en_nodos():
    """Acumula la masa de las barras en los nodos correspondientes."""
    global barras, masa_acumulada_por_nodo, A_Rope, A_Large, A_Medium, A_Small, gamma_fibra_carbono
    todos_los_nodos = set()
    for element in barras:
        _, nodo_i, nodo_j, _ = element
        todos_los_nodos.update([nodo_i, nodo_j])
    for nodo in todos_los_nodos:
        masa_nodo = 0
        for element in barras:
            _, nodo_i, nodo_j, Seccion = element
            if nodo in [nodo_i, nodo_j]:
                coord_i = np.array(ops.nodeCoord(nodo_i))
                coord_j = np.array(ops.nodeCoord(nodo_j))
                longitud = np.linalg.norm(coord_j - coord_i)
                area = {'S': A_Small, 'M': A_Medium, 'L': A_Large, 'Rope': A_Rope}[Seccion]
                masa_barra = longitud * area * gamma_fibra_carbono
                masa_nodo += masa_barra / 2
        masa_acumulada_por_nodo[nodo] = masa_nodo
    return masa_acumulada_por_nodo

def acumular_masa_paneles():
    """Acumula la masa de los paneles en los nodos correspondientes."""
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
    """Asigna las masas acumuladas a los nodos en OpenSees."""
    global masa_acumulada_por_nodo
    for nodo, masa in masa_acumulada_por_nodo.items():
        ops.mass(nodo, masa, masa, masa)

# Funciones de análisis
def realizar_analisis_frecuencias(num_modes):
    """Realiza un análisis modal y retorna las frecuencias naturales."""
    eigenvalues = ops.eigen(num_modes)
    eigenfrequencies = np.sqrt(eigenvalues) / (2 * np.pi)
    print(f'\nLa frecuencia natural más baja es: {eigenfrequencies[0]} Hz\n')
    return eigenfrequencies

def realizar_analisis_inercial(aceleraciones):
    """Realiza el análisis inercial para las aceleraciones dadas."""
    global barras
    fuerzas_maximas = {}
    solucion = 0

    for a in aceleraciones:
        solucion += 1
        # Definir la serie de tiempo y el patrón de carga
        ops.timeSeries('Constant', solucion)
        ops.pattern('Plain', solucion, solucion)

        # Aplicar las cargas nodales para todos los nodos
        for nodo, masa in masa_acumulada_por_nodo.items():
            F_inertial_x = masa * a[0]  # Fuerza inercial en X
            F_inertial_y = masa * a[1]  # Fuerza inercial en Y
            F_inertial_z = masa * a[2]  # Fuerza inercial en Z
            ops.load(nodo, F_inertial_x, F_inertial_y, F_inertial_z)

        # Configurar y ejecutar el análisis
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.algorithm('Linear')
        ops.analysis('Static')
        ops.analyze(1)

        # Obtener las fuerzas internas en las barras después del análisis
        for element in barras:
            element_id, nodo_i, nodo_j, _ = element
            # Obtener la fuerza interna en la barra
            fuerzas = ops.eleResponse(element_id, 'axialForce')
            fuerza_axial = fuerzas[0]  # Fuerza axial en el elemento
            # Acumular la suma de los cuadrados de las fuerzas axiales
            if element_id not in fuerzas_maximas:
                fuerzas_maximas[element_id] = fuerza_axial ** 2
            else:
                fuerzas_maximas[element_id] += fuerza_axial ** 2

        # Limpiar cargas y análisis para la siguiente iteración
        ops.remove('loadPattern', solucion)
        ops.wipeAnalysis()

    # Calcular la raíz cuadrada de la suma de los cuadrados (SRSS)
    for key in fuerzas_maximas:
        fuerzas_maximas[key] = (fuerzas_maximas[key]) ** 0.5

    return fuerzas_maximas

def realizar_analisis_termico():
    """Realiza el análisis térmico aplicando cargas térmicas a los elementos afectados."""
    global barras, barra_actual, alpha_T, deltaT

    thermal_load = alpha_T * deltaT
    # Crear un nuevo material que incluye la deformación térmica
    matTag_with_thermal = 3
    ops.uniaxialMaterial('InitStrainMaterial', matTag_with_thermal, 1, thermal_load)

    # Obtener los nodos de los paneles
    nodos_paneles_set = set()
    for panel in conexiones_paneles:
        for nodo in panel:
            nodos_paneles_set.add(nodo)

    # Filtrar las barras que conectan nodos en 'nodos_paneles'
    barras_paneles = []
    for barra in barras:
        barra_id, nodo_i, nodo_j, tamaño = barra
        if nodo_i in nodos_paneles_set and nodo_j in nodos_paneles_set:
            barras_paneles.append(barra)

    # Aplicar la carga térmica a las barras filtradas
    for barra in barras_paneles:
        # Reemplazar el elemento con uno que tenga el material con deformación térmica
        if barra[-1] == 'Rope':
            area = A_Rope
        elif barra[-1] == 'S':
            area = A_Small
        elif barra[-1] == 'M':
            area = A_Medium
        elif barra[-1] == 'L':
            area = A_Large

        # Reemplazar el elemento con el nuevo material
        ops.element('Truss', barra_actual, barra[1], barra[2], area, matTag_with_thermal, '-rho', gamma_fibra_carbono)
        barra_actual += 1

    # Configurar y ejecutar el análisis
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

# Funciones para cálculo de vectores normales y ángulos
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

# Funciones de visualización
def visualizar_panel_3_puntos(nodos_paneles, plotter, color='green'):
    """Visualiza paneles definidos por tres puntos."""
    for conexiones in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in conexiones])
        if len(puntos_panel) != 3:
            raise ValueError("Cada panel debe estar definido por exactamente 3 nodos.")
        surface = pv.PolyData(puntos_panel)
        surface.faces = np.hstack([[3], [0, 1, 2]])
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.3)

def visualizar_panel_4_puntos(nodos_paneles, plotter, color='green'):
    """Visualiza paneles definidos por cuatro puntos."""
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=1)

def visualizar_caja_satelite(caras_caja, barras, conexiones_paneles):
    """Visualiza la estructura completa de la caja satélite."""
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

def visualizar_desplazamiento_termico(barras, conexiones_paneles, escala):
    """Visualiza la estructura original y deformada térmicamente."""
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

def graficar_mapa_calor_inercia (barras, fuerzas_maximas):
    plotter = pv.Plotter()
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])

    f_max = max(fuerzas_maximas.values())

    lines = []
    colors = []

    fuerzas_normalizadas = {}

    for key in fuerzas_maximas:
        fuerzas_normalizadas[key] = fuerzas_maximas[key]/(f_max)

    # Crear lista de líneas y colores
    lines = []
    colors = []
    for element in barras:
        element_id, nodo_i, nodo_j, _ = element
        lines.append([2, nodo_i - 1, nodo_j - 1])
        
        fuerza_normalizada = fuerzas_normalizadas.get(element_id, 0)
        
        colors.append(255*fuerza_normalizada)

    truss = pv.PolyData(nodes)
    truss.lines = np.hstack(lines)

    for i, line in enumerate(truss.lines.reshape(-1, 3)):
        start_idx = line[1]
        end_idx = line[2]
        points = nodes[[start_idx, end_idx]]
        segment = pv.Line(points[0], points[1])
        
        # Determinar el color en la escala azul-verde-rojo
        color_value = colors[i]
        if color_value > 128:  # De verde a azul
            green_component = int(255 - 2 * (color_value - 128))  # Decrece de 255 a 0
            blue_component = int(2 * (color_value - 128))         # Aumenta de 0 a 255
            color_rgb = [0, green_component, blue_component]
        else:  # De rojo a verde
            red_component = int(255 - 2 * color_value)            # Decrece de 255 a 0
            green_component = int(2 * color_value)                # Aumenta de 0 a 255
            color_rgb = [red_component, green_component, 0]

        plotter.add_mesh(segment, color=color_rgb, line_width=2)

    # Agregar etiquetas de números de nodos (opcional)
    node_tags = ops.getNodeTags()
    labels = [str(tag) for tag in node_tags]
    plotter.add_point_labels(nodes, labels, point_size=5, font_size=12, text_color='black', name='labels')
    #visualizar_estructura_rigida(nodos_paneles, plotter, color='gold')
    #plotter.add_text(titulo, position='upper_left', font_size=10)
    plotter.add_legend([("Tension Maxima", "blue"), ("Tension Media", "green"),("Tension Minima", "red")])
    
    # Agregar la leyenda de la escala de colores en la parte inferior
    # Agregar la barra de colores con el rango de valores de 0 a f_max
    plotter.add_scalar_bar(title="Esfuerzo Interno Normalizado", 
                        position_x=0.3, position_y=0.05, 
                        width=0.4, height=0.05, 
                        label_font_size=10, title_font_size=12, 
                        n_labels=5)

    # Ajustar el rango de valores de la barra de colores
    plotter.update_scalar_bar_range([0, f_max*100])

    plotter.show_axes()
    plotter.show()



# Función principal
def main():
    global gamma_fibra_carbono, E_fibra_carbono
    global A_Rope, A_Small, A_Medium, A_Large
    global barras, conexiones_paneles, nodo_actual, barra_actual
    global alpha_T, deltaT

    # Inicialización y definición de materiales
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

    # Definir nodos de apoyo y generar barras
    nodes_apoyos_list = [
        np.array([[3.5, 3.3, 1], [0.3, 3.3, -1], [3.9, 0.9, -1]]),
        np.array([[-3.5, 3.3, 1], [-3.9, 0.9, -1], [-0.3, 3.3, -1]]),
        np.array([[-3.5, -3.3, 1], [-0.3, -3.3, -1], [-3.9, -0.9, -1]]),
        np.array([[3.5, -3.3, 1], [3.9, -0.9, -1], [0.3, -3.3, -1]])
    ]
    nodos_barras = [[] for _ in range(4)]
    for nodes_apoy, nodos_barra in zip(nodes_apoyos_list, nodos_barras):
        nodos_apoyos(nodes_apoy, nodos_barra)

    alphas = [45, 135, 225, 315]  # Ángulos para las barras
    for alpha, nodos_barra in zip(alphas, nodos_barras):
        nodos_base_barras([largo_inicial_barras, alpha, 0], nodos_barra, alpha)
        alargar_barras(alpha, nodos_barra)

    # Definir conexiones de paneles
    for i in range(4):
        nodos_paneles(nodos_barras[i], nodos_barras[(i + 1) % 4])

    # Paneles adicionales
    conexiones_paneles.append([nodos_barras[0][1][0], nodos_barras[1][1][0], nodos_barras[2][1][0]])
    conexiones_paneles.append([nodos_barras[2][1][0], nodos_barras[3][1][0], nodos_barras[0][1][0]])

    visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

    # Acumular masas y asignar masas a los nodos
    acumular_masa_barras_en_nodos()
    masa_barras = sum(masa_acumulada_por_nodo.values())
    area_paneles, masa_paneles = acumular_masa_paneles()
    print(f'Los paneles tienen un área total de {area_paneles} y una masa total de {masa_paneles}')
    if area_paneles < 3333.333333:
        raise ValueError("El área total de los paneles no es suficiente.")

    print(f'La masa total de las barras es: {masa_barras}')
    print(f'La proporción es: {masa_barras / masa_paneles}')
    asignar_masas_a_nodos()

    # Realizar análisis de frecuencias
    eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)

    # Análisis Inercial
    aceleracion_magnitud = 0.1 * 9.81
    aceleraciones = [
        [aceleracion_magnitud, 0, 0],
        [0, aceleracion_magnitud, 0],
        [0, 0, aceleracion_magnitud]
    ]
    fuerzas_maximas = realizar_analisis_inercial(aceleraciones)
    graficar_mapa_calor_inercia(barras, fuerzas_maximas)
    

    # Análisis Térmico
    realizar_analisis_termico()

    # Verificar ángulos entre vectores normales de los paneles
    lista_normal_inicial, lista_normal_final = calcular_vectores_normales_paneles(conexiones_paneles)
    for i in range(len(lista_normal_inicial)):
        angulo = calcular_angulo_entre_vectores(lista_normal_inicial[i], lista_normal_final[i])
        if angulo > 2:
            raise ValueError(f"El ángulo entre los vectores normales del panel es mayor a 2 grados")

    # Visualizar la estructura termica
    visualizar_desplazamiento_termico(barras, conexiones_paneles, escala=1000)

if __name__ == "__main__":
    main()
