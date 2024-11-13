import numpy as np
import openseespy.opensees as ops
import pyvista as pv
from math import atan

# Declaración de variables globales
nodo_actual, barra_actual = 0, 0
gamma_fibra_carbono = 1.91 * 1000
E_fibra_carbono = 338e9
gamma_panel = 1.1
D1_Rope, D2_Rope, A_Rope = 0.004, 0.000, None
D1_Small, D2_Small, A_Small = 0.005, 0.003, None
D1_Medium, D2_Medium, A_Medium = 0.1, 0.0995, None
D1_Large, D2_Large, A_Large = 0.15, 0.1485, None
ancho_barras, alto_barras, largo_inicial_barras = 2, 2, 5.5
largo_barras, espaciamiento = 30, 6
conexiones_paneles = []
masa_acumulada_por_nodo = {}
altura_torre, espaciado_torre = 6, 2
barras = []

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
    alpha_1 = 0
    alpha_2 = np.degrees(atan(1.65/ 1.95))
    alpha_3 = 90
    alpha_4 = alpha_2 + np.degrees(2 * atan(1.95/ 1.65))
    alpha_5 = 180
    alpha_6 = 2 * alpha_2 + alpha_4
    alpha_7 = 270
    alpha_8 = alpha_6 + np.degrees(2 * atan(1.95/ 1.65))
    return [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, alpha_7, alpha_8]

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
        if i == 0:
            ops.fix(nodo_actual, 1, 1, 1)

        else:
            ops.fix(nodo_actual, 1, 1, 0)
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
        conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'M')

def visualizar_panel_4_puntos(nodos_paneles, plotter, color='green'):
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=1)

def visualizar_panel_3_puntos(nodos_paneles, plotter, color='green'):
    for conexiones in nodos_paneles:
        # Obtener las coordenadas de los nodos del panel
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in conexiones])
        
        # Verificar que el panel tenga exactamente 3 puntos
        if len(puntos_panel) != 3:
            raise ValueError("Cada panel debe estar definido por exactamente 3 nodos.")
        
        # Crear el objeto PolyData para el panel triangular
        surface = pv.PolyData(puntos_panel)
        
        # Definir la cara del panel (triángulo) con 3 puntos
        surface.faces = np.hstack([[3], [0, 1, 2]])
        
        # Agregar el panel al plotter
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.3)

def conectar_barras(barra_1, barra_2, num_nodos_intermedios, nodos_barra_nueva, desfase_z=2, desfase_transversal=1):
        global nodo_actual

        
        print(barra_1)
        print(barra_2)
        # Obtener los nodos inicial y final de las barras
        nodo_inicial = barra_1[-1][1]
        nodo_final = barra_2[-1][2]
        
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


def visualizar_caja_satelite(caras_caja, barras, conexiones_paneles):
    # Obtener los tags y coordenadas de los nodos
    node_tags = ops.getNodeTags()
    node_coords = [ops.nodeCoord(tag) for tag in node_tags]
    nodes = np.array(node_coords)
    
    # Crear el plotter y agregar los nodos
    plotter = pv.Plotter()
    points = pv.PolyData(nodes)
    labels = {i + 1: str(tag) for i, tag in enumerate(node_tags)}
    plotter.add_point_labels(points, labels, font_size=10, text_color="red", always_visible=True)
    
    # Visualizar las caras de la caja y los paneles
    visualizar_panel_4_puntos(caras_caja, plotter, color='blue')
    visualizar_panel_3_puntos(conexiones_paneles, plotter, color='gold')
    
    # Definir el mapeo de colores para cada tamaño
    color_mapping = {
        'S': 'green',    # Barras pequeñas
        'M': 'orange',   # Barras medianas
        'L': 'red',      # Barras grandes
        'Rope': 'purple' # Otros tipos de barras
    }
    
    # Mantener registro de los tamaños ya agregados a la leyenda
    sizes_plotted = set()
    
    # Agregar las barras al plot con el color correspondiente
    for barra in barras:
        # Extraer información de la barra
        _, nodo_i, nodo_j, Tamaño = barra
        nodo_inicio, nodo_fin = nodo_i - 1, nodo_j - 1  # Ajustar índices para numpy
        
        # Obtener las coordenadas de los nodos
        coord_inicio = nodes[nodo_inicio]
        coord_fin = nodes[nodo_fin]
        linea = pv.Line(coord_inicio, coord_fin)
        
        # Obtener el color según el tamaño
        color_barra = color_mapping.get(Tamaño, 'black')
        
        # Agregar la barra al plotter
        if Tamaño not in sizes_plotted:
            # Agregar etiqueta para la leyenda
            plotter.add_mesh(linea, color=color_barra, line_width=2, label=Tamaño)
            sizes_plotted.add(Tamaño)
        else:
            # Agregar la barra sin etiqueta
            plotter.add_mesh(linea, color=color_barra, line_width=2)
    
    # Agregar la leyenda al plotter
    plotter.add_legend()
    
    # Mostrar ejes y renderizar
    plotter.show_axes()
    plotter.show()



def conectar_capas(nodos_capa_1, nodos_capa_2, Tamaño_Longitudinal, Tamaño_diagonal):
    global barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    area_longitudinal = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño_Longitudinal]
    area_diagonal = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño_diagonal]
    for i in range(len(nodos_capa_1)):
        if i == 0:
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], 'L'])
            barra_actual += 1
        else:
            ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], area_longitudinal, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i], Tamaño_Longitudinal])
            barra_actual += 1
    for i in range(len(nodos_capa_1)):
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], area_diagonal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], Tamaño_diagonal])
        barra_actual += 1
    for i in range(len(nodos_capa_1)):
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], area_diagonal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], Tamaño_diagonal])
        barra_actual += 1

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
                if Seccion == 'S':
                    area = A_Small
                elif Seccion == 'M':
                    area = A_Medium
                elif Seccion == 'L':
                    area = A_Large
                elif Seccion == 'Rope':
                    area = A_Rope
                    
                masa_barra = longitud * area * gamma_fibra_carbono
                masa_nodo += masa_barra / 2
        masa_acumulada_por_nodo[nodo] = masa_nodo
    return masa_acumulada_por_nodo

def area_tres_nodos_por_numero(nodo1, nodo2, nodo3):
    # Obtener las coordenadas de los nodos
    A = np.array(ops.nodeCoord(nodo1))
    B = np.array(ops.nodeCoord(nodo2))
    C = np.array(ops.nodeCoord(nodo3))
    
    # Calcular los vectores AB y AC
    AB = B - A
    AC = C - A
    
    # Calcular el área del triángulo ABC usando el producto cruzado
    area_ABC = 0.5 * np.linalg.norm(np.cross(AB, AC))
    
    return area_ABC

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
            if i == 1: #Estoy en el nodo de en medio
                masa_acumulada_por_nodo[nodos_conectados[i]] = masa_acumulada_por_nodo.get(nodos_conectados[i], 0) + (masa_por_panel*2)
            else:
                masa_acumulada_por_nodo[nodos_conectados[i]] = masa_acumulada_por_nodo.get(nodos_conectados[i], 0) + masa_por_panel
    return area_total, masa_total

def asignar_masas_a_nodos():
    global masa_acumulada_por_nodo
    for nodo, masa in masa_acumulada_por_nodo.items():
        ops.mass(nodo, masa, masa, masa)

def alargar_barras(alpha, nodos_barra):
    global nodo_actual
    for a in range(int(largo_barras/espaciamiento)):
        largo_adicional = (a+1)*espaciamiento
        x, y, z = coordenadas_cartesianas(largo_inicial_barras + largo_adicional, alpha, 0)
        x_1, y_1 = calcular_nuevo_punto_transversal(alpha, largo_inicial_barras + largo_adicional, ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(alpha, largo_inicial_barras + largo_adicional, -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z+alto_barras/2)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z - alto_barras/2)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z - alto_barras/2)
        nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])
        
        for i in range(len(nodos_barra[-1])):
            j = (i + 1) % len(nodos_barra[-1])
            generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'S')
        
        if a < (int(largo_barras/espaciamiento) /2) -1:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'M')

        else:
            conectar_capas(nodos_barra[-2], nodos_barra[-1], 'M', 'S')

def nodos_paneles(nodos_barra_A, nodos_barra_B):
    for i in range(len(nodos_barra_A)-2):
        j = i + 1
        #conexiones_paneles.append([nodos_barra_A[j][0], nodos_barra_B[j][0], nodos_barra_B[j+1][0], nodos_barra_A[j+1][0]])
        conexiones_paneles.append([nodos_barra_A[j][0], nodos_barra_B[j][0], nodos_barra_B[j+1][0]])
        conexiones_paneles.append([nodos_barra_B[j+1][0], nodos_barra_A[j+1][0], nodos_barra_A[j][0]])

def realizar_analisis_frecuencias(num_modes):
    eigenvalues = ops.eigen(num_modes)
    eigenfrequencies = np.sqrt((eigenvalues)) / (2 * np.pi)
    print('')
    print(f'La frecuencia natural más baja es: {eigenfrequencies[0]} Hz')
    print('')
    return eigenfrequencies

# Funciones relacionadas con la torre y funciones adicionales

def generar_torre(nodos_torre, up, Tamaño):
    global nodo_actual, barra_actual, A_Small, A_Medium, A_Large, gamma_fibra_carbono, altura_torre, espaciado_torre
    for i in range(int(altura_torre/espaciado_torre)):
        if up:
            a = 1
        else:
            a = -1
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
        if len(nodos_torre) > 1:
            conectar_capas(nodos_torre[-2], nodos_torre[-1], 'M', 'S')
        for i in range(len(nodos_torre)):
            if i != 0:
                for a in range(len(nodos_torre[i])):
                    j = (a+2)%len(nodos_torre[i])
                    ops.element('Truss', barra_actual, nodos_torre[i][a], nodos_torre[i][j], A_Rope, 2, '-rho', gamma_fibra_carbono)
                    barras.append([barra_actual, nodos_torre[i][a], nodos_torre[i][j], 'Rope'])
                    barra_actual += 1
    return nodos_torre

def conectar_mitades_torre_con_barras(torre, barras_grupo, conexion_final, Tamaño):
    global barra_actual, barras, A_Small, A_Medium, A_Large, A_Rope, gamma_fibra_carbono
    mitad_barra = int(round((len(barras_grupo[0][0]) - 1) / 2))
    mitad_torre = int(round((len(torre) - 1) / 2))
    if Tamaño == 'Rope':
        area = A_Rope
    elif Tamaño == 'S':
        area = A_Small
    elif Tamaño == 'M':
        area = A_Medium
    elif Tamaño == 'L':
        area = A_Large
    for esquina, barra_set in enumerate(barras_grupo):
        for barra in barra_set:
            nodo_mitad_barra = barra[mitad_barra][0]
            nodo_mitad_torre = torre[mitad_torre][esquina]
            #ops.element('Truss', barra_actual, nodo_mitad_barra, nodo_mitad_torre, area_rope, 1, '-rho', gamma_fibra_carbono)
            #barras.append([barra_actual, nodo_mitad_barra, nodo_mitad_torre, Tamaño])
            #barra_actual += 1
            nodo_final_barra = barra[conexion_final][0]
            nodo_final_torre = torre[-1][esquina]
            ops.element('Truss', barra_actual, nodo_final_barra, nodo_final_torre, area, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodo_final_barra, nodo_final_torre, Tamaño])
            barra_actual += 1

def main():
    global gamma_fibra_carbono, E_fibra_carbono, A_Rope, A_Small, A_Medium, A_Large
    global D1_Rope, D2_Rope, D1_Small, D2_Small, D1_Medium, D2_Medium, D1_Large, D2_Large
    global barras, conexiones_paneles, nodo_actual, barra_actual
    inicializar_modelo()
    definir_material(material_id=1, E=E_fibra_carbono)
    A_Rope = definir_seccion_tubo(D1=D1_Rope, D2=D2_Rope)
    A_Small = definir_seccion_tubo(D1=D1_Small, D2=D2_Small)
    A_Medium = definir_seccion_tubo(D1=D1_Medium, D2=D2_Medium)
    A_Large = definir_seccion_tubo(D1=D1_Large, D2=D2_Large)
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
    caras_caja = [
        [1,4,8,5],
        [2,3,7,6],
        [1,2,6,5],
        [3,4,8,7],
        [1,2,3,4],
        [5,6,7,8]
    ]
    nodes_apoyos_list = [
        np.array([[3.3, 0, 1], [3.9, 3, -1], [3.9, -3, -1]]),
        np.array([[3.9, 3.3, 1], [0.3, 3.3, -1], [3.9, 0.9, -1]]),
        np.array([[0, 3.3, 1], [-3, 3.3, -1], [3, 3.3, -1]]),
        np.array([[-3.9, 3.3, 1], [-3.9, 0.9, -1], [-0.3, 3.3, -1]]),
        np.array([[-3.9, 0, 1], [-3.9, -3, -1], [-3.9, 3, -1]]),
        np.array([[-3.9, -3.3, 1], [-0.3, -3.3, -1], [-3.9, -0.9, -1]]),
        np.array([[0, -3.3, 1], [3, -3.3, -1], [-3, -3.3, -1]]),
        np.array([[3.9, -3.3, 1], [3.9, -0.9, -1], [0.3, -3.3, -1]])
    ]
    nodos_barras = [[] for _ in range(8)]
    for nodes_apoy, nodos_barra in zip(nodes_apoyos_list, nodos_barras):
        nodos_apoyos(nodes_apoy, nodos_barra)
    alphas = angulos_barras()
    for alpha, nodos_barra in zip(alphas, nodos_barras):
        nodos_base_barras([largo_inicial_barras, alpha, 0], nodos_barra, alpha)
        alargar_barras(alpha, nodos_barra)
    for i in range(8):
        nodos_paneles(nodos_barras[i], nodos_barras[(i+1)%8])
    # Generación de la torre
    nodos_torre_1 = []
    nodos_base_torre_1 = np.array([
        [1, 0, 1.3],
        [0, 1, 1.3],
        [-1, 0, 1.3],
        [0, -1, 1.3]
    ])
    for coord in nodos_base_torre_1:
        nodo_actual += 1
        ops.node(nodo_actual, *coord)
        ops.fix(nodo_actual, 1, 1, 1)

    nodos_torre_1.append([nodo_actual - 3, nodo_actual - 2, nodo_actual - 1, nodo_actual])
    nodos_torre_1 = generar_torre(nodos_torre_1, up=True, Tamaño='L')

    nodos_torre_2 = []
    nodos_base_torre_2 = np.array([
        [1, 0, -1.3],
        [0, 1, -1.3],
        [-1, 0, -1.3],
        [0, -1, -1.3]
    ])
    for coord in nodos_base_torre_2:
        nodo_actual += 1
        ops.node(nodo_actual, *coord)
        ops.fix(nodo_actual, 1, 1, 1)



    nodos_torre_2.append([nodo_actual - 3, nodo_actual - 2, nodo_actual - 1, nodo_actual])
    nodos_torre_2 = generar_torre(nodos_torre_2, up=False, Tamaño='L')

    # Conexión de la torre con las barras
    barras_grupo_1 = [nodos_barras[7], nodos_barras[0], nodos_barras[1]]
    barras_grupo_2 = [nodos_barras[1], nodos_barras[2], nodos_barras[3]]
    barras_grupo_3 = [nodos_barras[3], nodos_barras[4], nodos_barras[5]]
    barras_grupo_4 = [nodos_barras[5], nodos_barras[6], nodos_barras[7]]
    conectar_mitades_torre_con_barras(nodos_torre_1, [barras_grupo_1, barras_grupo_2, barras_grupo_3, barras_grupo_4],-1, Tamaño='Rope')
    conectar_mitades_torre_con_barras(nodos_torre_2, [barras_grupo_1, barras_grupo_2, barras_grupo_3, barras_grupo_4],-3, Tamaño='Rope')
    #nodos_barra_1_2 = []
    #conectar_barras(nodos_barras[0], nodos_barras[1], 3, nodos_barra_1_2)

    #Que pasa si conecto nodos con cuerdas

    def conectar_brazos(barra_1, barra_2):
        global barra_actual, A_Medium, gamma_fibra_carbono
        nodo_1 = barra_1[-1][1]
        nodo_2 = barra_2[-4][0]

        ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Rope, 2, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodo_1, nodo_2, 'Rope'])
        barra_actual += 1

        nodo_1 = barra_1[-4][0]
        nodo_2 = barra_2[-1][2]

        ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Rope, 2, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodo_1, nodo_2, 'Rope'])
        barra_actual += 1
        

        

    conectar_brazos(nodos_barras[0], nodos_barras[1])
    conectar_brazos(nodos_barras[1], nodos_barras[2])
    conectar_brazos(nodos_barras[2], nodos_barras[3])
    conectar_brazos(nodos_barras[3], nodos_barras[4])
    conectar_brazos(nodos_barras[4], nodos_barras[5])
    conectar_brazos(nodos_barras[5], nodos_barras[6])
    conectar_brazos(nodos_barras[6], nodos_barras[7])
    conectar_brazos(nodos_barras[7], nodos_barras[0])





    # Acumulación de masas y análisis modal
    acumular_masa_barras_en_nodos()
    masa_barras = sum(masa_acumulada_por_nodo.values())
    area_paneles, masa_paneles = acumular_masa_paneles()
    print(f'Los paneles tienen un área total de {area_paneles} y una masa total de {masa_paneles}')
    if area_paneles < 3333.333333:
        raise ValueError("El área total de los paneles no es suficiente.")
        
    print(f'La masa total de las barras es: {masa_barras}')
    print(f'La proporción es: {(masa_barras/masa_paneles)}')
    asignar_masas_a_nodos()
    eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)
    visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

if __name__ == "__main__":
    main()
