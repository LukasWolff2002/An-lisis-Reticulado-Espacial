import numpy as np
import openseespy.opensees as ops
import pyvista as pv
from math import atan

#La torre me sube la frecuencia
#Que pasa si en vez de una barra, conecto como reticulado


# Declaración de variables globales
nodo_actual, barra_actual = 0, 0
gamma_fibra_carbono = 1.91 * 1000
E_fibra_carbono = 338e9
gamma_panel = 1.1
D1_Rope, D2_Rope, A_Rope = 0.01, 0.0005, None
D1_Small, D2_Small, A_Small = 0.05, 0.045, None
D1_Medium, D2_Medium, A_Medium = 0.008, 0.004, None
D1_Large, D2_Large, A_Large = 0.010, 0.006, None
ancho_barras, alto_barras, largo_inicial_barras = 0.5, 6, 6
largo_barras, espaciamiento = 120, 12
alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, alpha_7, alpha_8 = 0, 0, 0, 0, 0, 0,0, 0
nodos_barra_1, nodos_barra_2, nodos_barra_3, nodos_barra_4, nodos_barra_5, nodos_barra_6, nodos_barra_7, nodos_barra_8 = [], [], [], [], [], [], [], []
barras = []
conexiones_paneles = []
masa_acumulada_por_nodo = {}
nodos_torre_1, nodos_torre_2 = [], []
altura_torre, espaciado_torre = 20, 5
  

def inicializar_modelo():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

def definir_material(material_id, E):
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)

def definir_nodos(nodes):
    global nodo_actual
    for i, n in enumerate(nodes):
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))

def definir_nodos_fijos(nodes):
    global nodo_actual
    for i, n in enumerate(nodes):
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        ops.fix(nodo_actual, 1, 1, 1)
    

def definir_seccion_tubo(D1, D2):
    if D1 <= D2:
        raise ValueError("El diámetro exterior debe ser mayor que el diámetro interior.")
    A = np.pi * (D1**2 - D2**2) / 4
    return A

def angulos_barras():
    global alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, alpha_7, alpha_8
    alpha_1 = 0
    alpha_2 = np.degrees(atan(1.95 / 1.65))
    alpha_3 = 90
    alpha_4 = alpha_2 + np.degrees(2 * atan(1.65 / 1.95))
    alpha_5 = 180
    alpha_6 = 2 * alpha_2 + alpha_4
    alpha_7 = 270
    alpha_8 = alpha_6 + np.degrees(2 * atan(1.65 / 1.95))

def generar_elemento_axial(nodo_1, nodo_2, Tamaño):
    global barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    if Tamaño == 'S':
        ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Small, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodo_1, nodo_2, 'S'])
    elif Tamaño == 'M':
        ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Medium, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodo_1, nodo_2, 'M'])
    elif Tamaño == 'L':
        ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Large, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodo_1, nodo_2, 'L'])
    barra_actual += 1

def coordenadas_cartesianas(largo, angulo_xy, altura_z):
    # Convertir el ángulo de grados a radianes
    angulo_rad = np.radians(angulo_xy)

    # Calcular las coordenadas en x, y y z
    x = largo * np.cos(angulo_rad)
    y = largo * np.sin(angulo_rad)
    z = altura_z

    return x, y, z

def calcular_nuevo_punto_transversal(angulo_xy, largo, distancia_transversal):

    # Convertir el ángulo inicial a radianes
    angulo_rad = np.radians(angulo_xy)

    # Calcular las coordenadas del punto a lo largo de la línea generada por el ángulo inicial y el largo
    x_base = largo * np.cos(angulo_rad)
    y_base = largo * np.sin(angulo_rad)

    # Calcular el ángulo perpendicular para el desplazamiento transversal
    angulo_perpendicular_rad = angulo_rad + np.pi / 2  # Ángulo a 90 grados respecto a la línea inicial

    # Calcular las coordenadas del nuevo punto después del desplazamiento transversal
    x_nuevo = x_base + distancia_transversal * np.cos(angulo_perpendicular_rad)
    y_nuevo = y_base + distancia_transversal * np.sin(angulo_perpendicular_rad)

    return x_nuevo, y_nuevo

def nodos_apoyos (nodos_apoyo, nodos_barra):
    global nodo_actual
    
    for i, n in enumerate(nodos_apoyo):
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        #Dejo fijo el nodo
        if i == 0:
            ops.fix(nodo_actual, 1, 1, 1) #liberar z no cambia mucho el modelo
        else:
            ops.fix(nodo_actual, 1, 1, 1)

    nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

def nodos_base_barras(nodos_iniciales_barra, nodos_barra, alpha):
    global nodo_actual, barra_actual, ancho_barras, alto_barras, largo_inicial_barras
    nodos_iniciales_barra = [largo_inicial_barras, alpha, alto_barras/2]
    x, y, z = coordenadas_cartesianas(nodos_iniciales_barra[0], nodos_iniciales_barra[1], nodos_iniciales_barra[2])
    x_1, y_1 = calcular_nuevo_punto_transversal(nodos_iniciales_barra[1], nodos_iniciales_barra[0], ancho_barras)
    x_2, y_2 = calcular_nuevo_punto_transversal(nodos_iniciales_barra[1], nodos_iniciales_barra[0], -ancho_barras)
    
    nodo_actual += 1
    ops.node(nodo_actual, x, y, z)
    nodo_actual += 1
    ops.node(nodo_actual, x_1, y_1, z - alto_barras)
    nodo_actual += 1
    ops.node(nodo_actual, x_2, y_2, z - alto_barras)

    nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

    for i in range(len(nodos_barra[-1])):
        j = (i + 1) % len(nodos_barra[-1])
        #Las conexiones entre 3 nodos van a ser M
        generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'M')

    #Ahora conecto con los apoyos
    if len(nodos_barra) > 1:
        conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'S')

def visualizar_panel(nodos_paneles, plotter, color='green'):
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.7)

def visualizar_caja_satelite(caras_caja, barras, conexiones_paneles):
    # Obtener los tags de los nodos en el modelo de OpenSees
    node_tags = ops.getNodeTags()

    # Crear una lista para almacenar las coordenadas de cada nodo
    node_coords = []

    # Iterar sobre cada nodo y obtener sus coordenadas
    for tag in node_tags:
        coords = ops.nodeCoord(tag)  # Obtener las coordenadas (x, y, z) del nodo
        node_coords.append(coords)

    # Convertir la lista a un array de numpy
    nodes = np.array(node_coords)

    # Crear el objeto de visualización en PyVista
    plotter = pv.Plotter()
    points = pv.PolyData(nodes)
    
    # Agregar etiquetas con el número de nodo en cada punto y hacerlas siempre visibles
    labels = {i + 1: str(tag) for i, tag in enumerate(node_tags)}
    plotter.add_point_labels(points, labels, font_size=10, text_color="red", always_visible=True)

    # Visualizar las caras de la caja satélite como paneles
    visualizar_panel(caras_caja, plotter, color='blue')
    visualizar_panel(conexiones_paneles, plotter, color='gold')

    # Dibujar las barras que conectan los nodos
    for barra in barras:
        _, nodo_i, nodo_j, __ = barra  # Ignorar el número de barra, solo usar nodo_i y nodo_j
        nodo_inicio, nodo_fin = nodo_i - 1, nodo_j - 1  # Restar 1 para índices de numpy
        coord_inicio = nodes[nodo_inicio]
        coord_fin = nodes[nodo_fin]
        linea = pv.Line(coord_inicio, coord_fin)
        plotter.add_mesh(linea, color="black", line_width=2)  # Agregar barra entre los nodos

    # Mostrar ejes y renderizar la visualización
    plotter.show_axes()
    plotter.show()

def conectar_capas (nodos_capa_1, nodos_capa_2, Tamaño_Longitudinal, Tamaño_diagonal):
    global barra_actual, A_Small, gamma_fibra_carbono
    #Hago la conexión entre las capas
    for i in range(len(nodos_capa_1)):
        j = (i + 1) % len(nodos_capa_1)
        if Tamaño_Longitudinal == 'S':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], A_Small, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], 'S'])
        elif Tamaño_Longitudinal == 'M':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], 'M'])
        elif Tamaño_Longitudinal == 'L':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], 'L'])
        barra_actual += 1

    #Agrego la conexion diagonal
    for i in range(len(nodos_capa_1)):
        if Tamaño_diagonal == 'S':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], A_Small, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], 'S'])
        elif Tamaño_diagonal == 'M':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], 'M'])
        elif Tamaño_diagonal == 'L':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], 'L'])
        barra_actual += 1
    #Agrego la conexion diagonal en el otro sentido
    for i in range(len(nodos_capa_1)):
        if Tamaño_diagonal == 'S':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], A_Small, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], 'S'])
        elif Tamaño_diagonal == 'M':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], 'M'])
        elif Tamaño_diagonal == 'L':
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], 'L'])
        barra_actual += 1



def acumular_masa_barras_en_nodos():
    global barras, masa_acumulada_por_nodo

    #Actualmente estoy asumiendo el area minima para cada barra

    # Crear un diccionario para almacenar la masa acumulada por nodo
    

    # Obtener la lista de todos los nodos en el modelo
    todos_los_nodos = set()
    for element in barras:
        _, nodo_i, nodo_j, Seccion = element
        todos_los_nodos.add(nodo_i)
        todos_los_nodos.add(nodo_j)

    # Para cada nodo, encontrar las barras conectadas a él
    for nodo in todos_los_nodos:
        masa_nodo = 0
        # Buscar elementos conectados al nodo actual
        for element in barras:
            element_id, nodo_i, nodo_j, Seccion = element
            if nodo == nodo_i or nodo == nodo_j:
                # Calcular la masa de la barra
                coord_i = np.array(ops.nodeCoord(nodo_i))
                coord_j = np.array(ops.nodeCoord(nodo_j))
                longitud = np.linalg.norm(coord_j - coord_i)
                if Seccion == 'S':
                    masa_barra = longitud * A_Small * gamma_fibra_carbono
                elif Seccion == 'M':
                    masa_barra = longitud * A_Medium * gamma_fibra_carbono
                elif Seccion == 'L':
                    masa_barra = longitud * A_Large * gamma_fibra_carbono
                
                # Agregar la mitad de la masa de la barra al nodo
                masa_nodo += masa_barra / 2
        # Almacenar la masa acumulada en el nodo
        masa_acumulada_por_nodo[nodo] = masa_nodo

    return masa_acumulada_por_nodo

def area_cuatro_nodos_por_numero(nodo1, nodo2, nodo3, nodo4):
    
    # Obtener las coordenadas de cada nodo
    A = np.array(ops.nodeCoord(nodo1))
    B = np.array(ops.nodeCoord(nodo2))
    C = np.array(ops.nodeCoord(nodo3))
    D = np.array(ops.nodeCoord(nodo4))
    
    # Calcular el área del triángulo ABC usando el producto cruzado
    AB = B - A
    AC = C - A
    area_ABC = 0.5 * np.linalg.norm(np.cross(AB, AC))
    
    # Calcular el área del triángulo ACD usando el producto cruzado
    AD = D - A
    area_ACD = 0.5 * np.linalg.norm(np.cross(AC, AD))
    
    # Sumamos las áreas de los dos triángulos para obtener el área total
    area_total = area_ABC + area_ACD
    return area_total

def acumular_masa_paneles():
    area_total = 0
    masa_total = 0
    global masa_acumulada_por_nodo, conexiones_paneles, gamma_panel
    for nodos_conectados in conexiones_paneles:
        area_panel = area_cuatro_nodos_por_numero(*nodos_conectados)
        area_total += area_panel
        masa_panel = area_panel * gamma_panel
        masa_total += masa_panel
        masa_por_panel = masa_panel / 4
        for nodo in nodos_conectados:
            #print(f"Nodo {nodo} - Masa asignada (Barras): {masa_acumulada_por_nodo.get(nodo, 0)}")
            masa_acumulada_por_nodo[nodo] = masa_acumulada_por_nodo.get(nodo, 0) + masa_por_panel

    return area_total, masa_total

def asignar_masas_a_nodos():
    global masa_acumulada_por_nodo
    for nodo, masa in masa_acumulada_por_nodo.items():
        ops.mass(nodo, masa, masa, masa)

def alargar_barras (alpha, nodos_barra):
    global nodo_actual
    for i in range(int(largo_barras/espaciamiento)):
        largo_adicional = (i+1)*espaciamiento
        nodos_nuevos = [largo_inicial_barras + largo_adicional, alpha, 0]
        x, y, z = coordenadas_cartesianas(nodos_nuevos[0], nodos_nuevos[1], nodos_nuevos[2])
        x_1, y_1 = calcular_nuevo_punto_transversal(nodos_nuevos[1], nodos_nuevos[0], ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(nodos_nuevos[1], nodos_nuevos[0], -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z+alto_barras/2)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z - alto_barras/2)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z - alto_barras/2)

        nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

        for i in range(len(nodos_barra[-1])):
            j = (i + 1) % len(nodos_barra[-1])
            generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'M')

        conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'S')

def nodos_paneles(nodos_barra_A, nodos_barra_B):
    for i in range(len(nodos_barra_A)-2):
        j = i + 1
        conexiones_paneles.append([nodos_barra_A[j][1], nodos_barra_B[j][2], nodos_barra_B[j+1][2], nodos_barra_A[j+1][1]])

def main():
    global gamma_fibra_carbono, E_fibra_carbono, A_Rope, A_Small, A_Medium, A_Large
    global D1_Rope, D1_Rope, D1_Small, D2_Small, D1_Medium, D2_Medium, D1_Large, D2_Large
    global nodos_barra_1, nodos_barra_2, nodos_barra_3, nodos_barra_4, barras
    global largo_barras, espaciamiento, conexiones_paneles, nodo_actual, barra_actual
    global altura_torre, espaciado_torre, nodos_torre_1, nodos_torre_2

    inicializar_modelo()

    definir_material(material_id=1, E=E_fibra_carbono)

    A_Rope = definir_seccion_tubo(D1=D1_Rope, D2=D2_Rope)
    A_Small = definir_seccion_tubo( D1=D1_Small, D2=D2_Small)
    A_Medium = definir_seccion_tubo(D1=D1_Medium, D2=D2_Medium)
    A_Large = definir_seccion_tubo(D1=D1_Large, D2=D2_Large)


    nodes = np.array([
        [3.3, 3.9, 1.3],
        [3.3, -3.9, 1.3],
        [3.3, -3.9, -1.3],
        [3.3, 3.9, -1.3],
        [-3.3, 3.9, 1.3],
        [-3.3, -3.9, 1.3],
        [-3.3, -3.9, -1.3],
        [-3.3, 3.9, -1.3]
    ])
    definir_nodos_fijos(nodes)

    #Defino los nodos que definiran las caras de la caja satelite
    caras_caja = [
        [1,4,8,5],
        [2,3,7,6],
        [1,2,6,5],
        [3,4,8,7],
        [1,2,3,4],
        [5,6,7,8]
    ]

    nodes_apoy_1 = np.array([
        [3.3, 0, 1],
        [3.3, 2, -1],
        [3.3, -2, -1]])

    nodes_apoy_2 = np.array([
                [3.3, 3.9, 1],
                [0.3, 3.9, -1],
                [3.3, 0.9, -1]])

    nodes_apoy_3 = np.array([
                [0, 3.9, 1],
                [-2, 3.9, -1],
                [2, 3.9, -1]])

    nodes_apoy_4 = np.array([
                [-3.3, 3.9, 1],
                [-3.3, 0.9, -1],
                [-0.3, 3.9, -1]])

    nodes_apoy_5 = np.array([
                [-3.3, 0, 1],
                [-3.3, -2, -1],
                [-3.3, 2, -1]])

    nodes_apoy_6 = np.array([
                [-3.3, -3.9, 1],
                [-0.3, -3.9, -1],
                [-3.3, -0.9, -1]])
    
    nodes_apoy_7 = np.array([
                [0, -3.9, 1],
                [2, -3.9, -1],
                [-2, -3.9, -1]])

    nodes_apoy_8 = np.array([
                [3.3, -3.9, 1],
                [3.3, -0.9, -1],
                [0.3, -3.9, -1]])
    
    nodos_apoyos(nodes_apoy_1, nodos_barra_1)
    nodos_apoyos(nodes_apoy_2, nodos_barra_2)
    nodos_apoyos(nodes_apoy_3, nodos_barra_3)
    nodos_apoyos(nodes_apoy_4, nodos_barra_4)
    nodos_apoyos(nodes_apoy_5, nodos_barra_5)
    nodos_apoyos(nodes_apoy_6, nodos_barra_6)
    nodos_apoyos(nodes_apoy_7, nodos_barra_7)
    nodos_apoyos(nodes_apoy_8, nodos_barra_8)

    angulos_barras()

    nodos_base_barras([largo_inicial_barras, alpha_1, 0], nodos_barra_1, alpha_1)
    nodos_base_barras([largo_inicial_barras, alpha_2, 0], nodos_barra_2, alpha_2)
    nodos_base_barras([largo_inicial_barras, alpha_3, 0], nodos_barra_3, alpha_3)
    nodos_base_barras([largo_inicial_barras, alpha_4, 0], nodos_barra_4, alpha_4)
    nodos_base_barras([largo_inicial_barras, alpha_5, 0], nodos_barra_5, alpha_5)
    nodos_base_barras([largo_inicial_barras, alpha_6, 0], nodos_barra_6, alpha_6)
    nodos_base_barras([largo_inicial_barras, alpha_7, 0], nodos_barra_7, alpha_7)
    nodos_base_barras([largo_inicial_barras, alpha_8, 0], nodos_barra_8, alpha_8)

    alargar_barras(alpha_1, nodos_barra_1)
    alargar_barras(alpha_2, nodos_barra_2)
    alargar_barras(alpha_3, nodos_barra_3)
    alargar_barras(alpha_4, nodos_barra_4)
    alargar_barras(alpha_5, nodos_barra_5)
    alargar_barras(alpha_6, nodos_barra_6)
    alargar_barras(alpha_7, nodos_barra_7)
    alargar_barras(alpha_8, nodos_barra_8)


    nodos_paneles(nodos_barra_1, nodos_barra_2)
    nodos_paneles(nodos_barra_2, nodos_barra_3)
    nodos_paneles(nodos_barra_3, nodos_barra_4)
    nodos_paneles(nodos_barra_4, nodos_barra_5)
    nodos_paneles(nodos_barra_5, nodos_barra_6)
    nodos_paneles(nodos_barra_6, nodos_barra_7)
    nodos_paneles(nodos_barra_7, nodos_barra_8)
    nodos_paneles(nodos_barra_8, nodos_barra_1)


    #--------------------------------------------
    #Genero la torre 
    #--------------------------------------------

    nodos_base_torre_1 = np.array([
        [1, 0, 1.3],
        [0, 1, 1.3],
        [-1, 0, 1.3],
        [0, -1, 1.3]
    ])

    for i in range(len(nodos_base_torre_1)):
        x, y, z = nodos_base_torre_1[i]
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        ops.fix(nodo_actual, 1, 1, 1)

    nodos_torre_1 = [[nodo_actual - 3, nodo_actual - 2, nodo_actual - 1, nodo_actual]]

    nodos_base_torre_2 = np.array([
        [1, 0, -1.3],
        [0, 1, -1.3],
        [-1, 0, -1.3],
        [0, -1, -1.3]
    ])

    for i in range(len(nodos_base_torre_2)):
        x, y, z = nodos_base_torre_2[i]
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        ops.fix(nodo_actual, 1, 1, 1)
        
    
    nodos_torre_2 = [[nodo_actual - 3, nodo_actual - 2, nodo_actual - 1, nodo_actual]]



   

    def generar_torre(nodos_torre, up, Tamaño):
        global nodo_actual, barra_actual, A_Small, gamma_fibra_carbono, altura_torre, espaciado_torre


        for i in range(int(altura_torre/espaciado_torre)):
            #gnero siguiente capa de altura
            print('------------------')
            print(i)
        
            if up:
                a = 1
        
            else:
                a =-1
                
            altura = a*(i+1)*espaciado_torre
        
            nodo_actual += 1
            ops.node(nodo_actual, 1, 0, altura)
            nodo_actual += 1
            ops.node(nodo_actual, 0, 1, altura)
            nodo_actual += 1
            ops.node(nodo_actual, -1, 0, altura)
            nodo_actual += 1
            ops.node(nodo_actual, 0, -1, altura)

            nodos_torre.append([nodo_actual - 3, nodo_actual - 2, nodo_actual - 1, nodo_actual])

            for b in range(len(nodos_torre[-1])):
                j = (b + 1) % len(nodos_torre[-1])
                generar_elemento_axial(nodos_torre[-1][b], nodos_torre[-1][j], 'M')

            #Ahora conecto las capas
            conectar_capas(nodos_torre[-2], nodos_torre[-1], 'L', 'S')

        #Ahora agrego las conexiones x en las capas
        for i in range(len(nodos_torre)):
            if i != 0:
                if Tamaño == 'S':
                    ops.element('Truss', barra_actual, nodos_torre[i][0], nodos_torre[i][0], A_Small, 1, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][0], nodos_torre[i][0], 'S'])
                    barra_actual += 1
                    ops.element('Truss', barra_actual, nodos_torre[i][1], nodos_torre[i][1], A_Small, 1, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][1], nodos_torre[i][1], 'S'])
                    barra_actual += 1
                
                elif Tamaño == 'M':
                    ops.element('Truss', barra_actual, nodos_torre[i][0], nodos_torre[i][0], A_Medium, 1, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][0], nodos_torre[i][0], 'M'])
                    barra_actual += 1
                    ops.element('Truss', barra_actual, nodos_torre[i][1], nodos_torre[i][1], A_Medium, 1, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][1], nodos_torre[i][1], 'M'])
                    barra_actual += 1

                elif Tamaño == 'L':
                    ops.element('Truss', barra_actual, nodos_torre[i][0], nodos_torre[i][0], A_Large, 1, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][0], nodos_torre[i][0], 'L'])
                    barra_actual += 1
                    ops.element('Truss', barra_actual, nodos_torre[i][1], nodos_torre[i][1], A_Large, 1, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][1], nodos_torre[i][1], 'L'])
                    barra_actual += 1
          

        return nodos_torre
        
    nodos_torre_1 = generar_torre(nodos_torre_1, True, 'L')
    #nodos_torre_2 = generar_torre(nodos_torre_2, False)


    

    #Bien, ahora debo conetar la torre a los brazos

    
    
    def conectar_mitades_torre_con_barras(torre, barras_grupo, gamma_fibra_carbono, Tamaño):

        global barra_actual, barras, A_Rope

        # La mitad de la barra sin contar el espacio inicial
        mitad_barra = int(round((len(barras_grupo[0][0])) / 2, 0)) 
        print('la mitad de la barra es:', mitad_barra)
        mitad_torre = int(round((len(torre) - 1) / 2, 0)) 

        

        for esquina, barra_set in enumerate(barras_grupo):
            print(mitad_barra)
            for barra in barra_set:

                if Tamaño == 'S':
                    A_Rope = A_Small
                
                elif Tamaño == 'M':
                    A_Rope = A_Medium
                
                elif Tamaño == 'L':
                    A_Rope = A_Large

                # Conectar nodo de mitad de la barra con el nodo de mitad de la torre
                ops.element('Truss', barra_actual, barra[mitad_barra][0], torre[mitad_torre][esquina], A_Rope, 1, '-rho', gamma_fibra_carbono)
                barras.append([barra_actual, barra[mitad_barra][0], torre[mitad_torre][esquina], Tamaño])
                barra_actual += 1

                # Conectar nodo final de la barra con el nodo final de la torre
                ops.element('Truss', barra_actual, barra[-1][0], torre[-1][esquina], A_Rope, 1, '-rho', gamma_fibra_carbono)
                barras.append([barra_actual, barra[-1][0], torre[-1][esquina], Tamaño])
                barra_actual += 1

    # Ejemplo de uso:
    barras_grupo_1 = [nodos_barra_8, nodos_barra_1, nodos_barra_2]  # Barras para la primera esquina
    barras_grupo_2 = [nodos_barra_2, nodos_barra_3, nodos_barra_4]  # Barras para la segunda esquina
    barras_grupo_3 = [nodos_barra_4, nodos_barra_5, nodos_barra_6]  # Barras para la tercera esquina
    barras_grupo_4 = [nodos_barra_6, nodos_barra_7, nodos_barra_8]  # Barras para la cuarta esquina

    # Conectar cada esquina con sus barras correspondientes
    conectar_mitades_torre_con_barras(nodos_torre_1, [barras_grupo_1, barras_grupo_2, barras_grupo_3, barras_grupo_4], gamma_fibra_carbono, 'L')
    #conectar_mitades_torre_con_barras(nodos_torre_2, [barras_grupo_1, barras_grupo_2, barras_grupo_3, barras_grupo_4], A_Small, gamma_fibra_carbono)





    def conectar_barras(barra_1, nodo_1, barra_2,nodo_2, num_nodos_intermedios, nodos_barra_nueva, desfase_z=2, desfase_transversal=1):
        global nodo_actual

        # Obtener los nodos inicial y final de las barras
        #nodo_inicial = barra_1[3][1]
        nodo_inicial = barra_1[nodo_1][0]
        nodo_final = barra_2[-1][nodo_2]
        #Antes nodo 1 era 2
        
        # Obtener coordenadas de los nodos inicial y final
        coord_inicial = np.array(ops.nodeCoord(nodo_inicial))
        coord_final = np.array(ops.nodeCoord(nodo_final))
        
        # Calcular el incremento en cada dirección
        delta = (coord_final - coord_inicial) / (num_nodos_intermedios + 1)
        
        # Calcular el vector perpendicular para el desfase transversal
        direccion = coord_final - coord_inicial
        perpendicular = np.array([-direccion[1], direccion[0], 0])  # Perpendicular en el plano XY
        perpendicular = perpendicular / np.linalg.norm(perpendicular)  # Normalizar
        
        # Crear nodos intermedios
        for i in range(1, num_nodos_intermedios + 1):
            # Coordenadas del nodo intermedio sin desfases
            coord_nodo = coord_inicial + i * delta
            
            # Aplicar el desfase en z y el desfase transversal
            coord_nodo_desfasado_z = coord_nodo + np.array([0, 0, desfase_z]) 
            coord_nodo_desfazado_lat = coord_nodo + perpendicular * desfase_transversal
            cooed_nodo_desfasado_lat_2 = coord_nodo - perpendicular * desfase_transversal
            
    
            # Crear el nodo en OpenSees
            nodo_actual += 1
            ops.node(nodo_actual, *coord_nodo_desfasado_z)
            nodo_actual += 1
            ops.node(nodo_actual, *coord_nodo_desfazado_lat)
            nodo_actual += 1
            ops.node(nodo_actual, *cooed_nodo_desfasado_lat_2)

            nodos_barra_nueva.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

            #Conecto el nodo inicio con los primeros nodos
            if i == 1:
                for a in range(len(nodos_barra_nueva[-1])):
                    generar_elemento_axial(nodo_inicial, nodos_barra_nueva[-1][a], 'M')

            for a in range(len(nodos_barra_nueva[-1])):
     
                j = (a + 1) % len(nodos_barra_nueva[-1])
                generar_elemento_axial(nodos_barra_nueva[-1][a], nodos_barra_nueva[-1][j], 'M')

            if len(nodos_barra_nueva) > 1:
               
                conectar_capas(nodos_barra_nueva[-1], nodos_barra_nueva[-2], 'M', 'S')

        # Conectar el último nodo intermedio con el nodo final
        for a in range(len(nodos_barra_nueva[-1])):
            generar_elemento_axial(nodos_barra_nueva[-1][a], nodo_final, 'M')

            


    nodos_barra_1_2 = []
    nodos_barra_2_3 = []
    nodos_barra_3_4 = []
    nodos_barra_4_5 = []
    nodos_barra_5_6 = []
    nodos_barra_6_7 = []
    nodos_barra_7_8 = []
    nodos_barra_8_1 = []
    
    #conectar_barras(nodos_barra_5, -1, nodos_torre_1, 2, 10, nodos_barra_1_2)
    #conectar_barras(nodos_barra_7, -1, nodos_torre_1, 3, 10, nodos_barra_2_3)
    #conectar_barras(nodos_barra_1, -1, nodos_torre_1, 0, 10, nodos_barra_3_4)
    #conectar_barras(nodos_barra_3, -1, nodos_torre_1, 1, 10, nodos_barra_4_5)

    #conectar_barras(nodos_barra_1, nodos_barra_2, 3, nodos_barra_1_2)
    #conectar_barras(nodos_barra_2, nodos_barra_3, 3, nodos_barra_2_3)
    #conectar_barras(nodos_barra_3, nodos_barra_4, 3, nodos_barra_3_4)
    #conectar_barras(nodos_barra_4, nodos_barra_5, 3, nodos_barra_4_5)
    #conectar_barras(nodos_barra_5, nodos_barra_6, 3, nodos_barra_5_6)
    #conectar_barras(nodos_barra_6, nodos_barra_7, 3, nodos_barra_6_7)
    #conectar_barras(nodos_barra_7, nodos_barra_8, 3, nodos_barra_7_8)
    #conectar_barras(nodos_barra_8, nodos_barra_1, 3, nodos_barra_8_1)



    







        






    #visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

    


    #--------------------------------------------
    #Ahora agrego las masas a los nodos, para hacer los analisis modales
    #--------------------------------------------

    #En este punto ya tengo las masas de las barras, asumiendo un area S
    acumular_masa_barras_en_nodos()
    masa_barras = sum(masa_acumulada_por_nodo.values())
    area_paneles, masa_paneles = acumular_masa_paneles()
    print(f'los paneles tienen un area total de {area_paneles} y una masa total de {masa_paneles}')
    if area_paneles < 3333.33333:
        print('El area de los paneles es menor a 3333333.333')

    print(f'La masa total de las barras es: {masa_barras}')
    print(f'La proporcion es: {(masa_barras/masa_paneles)}')
    asignar_masas_a_nodos()

    #--------------------------------------------
    #Analisis de Frecuencias
    #--------------------------------------------

    def realizar_analisis_frecuencias(num_modes):
        eigenvalues = ops.eigen(num_modes)
        eigenfrequencies = np.sqrt((eigenvalues)) / (2 * np.pi)
        print(f'La frecuencia natural mas baja es: {eigenfrequencies} Hz')
        return eigenfrequencies

    eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)
    

   

    #ops.printModel('-node')

    visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

if __name__ == "__main__":
    main()
