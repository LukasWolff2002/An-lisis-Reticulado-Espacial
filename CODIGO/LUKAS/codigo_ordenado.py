import numpy as np
import openseespy.opensees as ops
import pyvista as pv
from math import atan

# Declaración de variables globales
nodo_actual, barra_actual = 0, 0
gamma_fibra_carbono = 1.91 * 1000
E_fibra_carbono = 338e9
gamma_panel = 1.1
D1_Rope, D2_Rope, A_Rope = 0.01, 0.0005, None
D1_Small, D2_Small, A_Small = 0.05, 0.045, None
D1_Medium, D2_Medium, A_Medium = 0.008, 0.004, None
D1_Large, D2_Large, A_Large = 0.010, 0.006, None
ancho_barras, alto_barras, largo_inicial_barras = 1, 2, 6
largo_barras, espaciamiento = 120, 12
conexiones_paneles = []
masa_acumulada_por_nodo = {}
altura_torre, espaciado_torre = 20, 5
barras = []

def inicializar_modelo():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

def definir_material(material_id, E):
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)

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
    alpha_2 = np.degrees(atan(1.95 / 1.65))
    alpha_3 = 90
    alpha_4 = alpha_2 + np.degrees(2 * atan(1.65 / 1.95))
    alpha_5 = 180
    alpha_6 = 2 * alpha_2 + alpha_4
    alpha_7 = 270
    alpha_8 = alpha_6 + np.degrees(2 * atan(1.65 / 1.95))
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
        generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'M')
    if len(nodos_barra) > 1:
        conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'S')

def visualizar_panel(nodos_paneles, plotter, color='green'):
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.7)

def visualizar_caja_satelite(caras_caja, barras, conexiones_paneles):
    node_tags = ops.getNodeTags()
    node_coords = [ops.nodeCoord(tag) for tag in node_tags]
    nodes = np.array(node_coords)
    plotter = pv.Plotter()
    points = pv.PolyData(nodes)
    labels = {i + 1: str(tag) for i, tag in enumerate(node_tags)}
    plotter.add_point_labels(points, labels, font_size=10, text_color="red", always_visible=True)
    visualizar_panel(caras_caja, plotter, color='blue')
    visualizar_panel(conexiones_paneles, plotter, color='gold')
    for barra in barras:
        _, nodo_i, nodo_j, __ = barra
        nodo_inicio, nodo_fin = nodo_i - 1, nodo_j - 1
        coord_inicio = nodes[nodo_inicio]
        coord_fin = nodes[nodo_fin]
        linea = pv.Line(coord_inicio, coord_fin)
        plotter.add_mesh(linea, color="black", line_width=2)
    plotter.show_axes()
    plotter.show()

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
    for i in range(len(nodos_capa_1)):
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], area_diagonal, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i-1)%len(nodos_capa_2)], Tamaño_diagonal])
        barra_actual += 1

def acumular_masa_barras_en_nodos():
    global barras, masa_acumulada_por_nodo
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
                area = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Seccion]
                masa_barra = longitud * area * gamma_fibra_carbono
                masa_nodo += masa_barra / 2
        masa_acumulada_por_nodo[nodo] = masa_nodo
    return masa_acumulada_por_nodo

def area_cuatro_nodos_por_numero(nodo1, nodo2, nodo3, nodo4):
    A = np.array(ops.nodeCoord(nodo1))
    B = np.array(ops.nodeCoord(nodo2))
    C = np.array(ops.nodeCoord(nodo3))
    D = np.array(ops.nodeCoord(nodo4))
    AB = B - A
    AC = C - A
    area_ABC = 0.5 * np.linalg.norm(np.cross(AB, AC))
    AD = D - A
    area_ACD = 0.5 * np.linalg.norm(np.cross(AC, AD))
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
            masa_acumulada_por_nodo[nodo] = masa_acumulada_por_nodo.get(nodo, 0) + masa_por_panel
    return area_total, masa_total

def asignar_masas_a_nodos():
    global masa_acumulada_por_nodo
    for nodo, masa in masa_acumulada_por_nodo.items():
        ops.mass(nodo, masa, masa, masa)

def alargar_barras(alpha, nodos_barra):
    global nodo_actual
    for i in range(int(largo_barras/espaciamiento)):
        largo_adicional = (i+1)*espaciamiento
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
            generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j], 'M')
        conectar_capas(nodos_barra[-2], nodos_barra[-1], 'L', 'S')

def nodos_paneles(nodos_barra_A, nodos_barra_B):
    for i in range(len(nodos_barra_A)-2):
        j = i + 1
        conexiones_paneles.append([nodos_barra_A[j][1], nodos_barra_B[j][2], nodos_barra_B[j+1][2], nodos_barra_A[j+1][1]])

def realizar_analisis_frecuencias(num_modes):
    eigenvalues = ops.eigen(num_modes)
    eigenfrequencies = np.sqrt((eigenvalues)) / (2 * np.pi)
    print(f'La frecuencia natural más baja es: {eigenfrequencies[0]} Hz')
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
            conectar_capas(nodos_torre[-2], nodos_torre[-1], 'L', 'S')
    return nodos_torre

def conectar_mitades_torre_con_barras(torre, barras_grupo, Tamaño):
    global barra_actual, barras, A_Small, A_Medium, A_Large, gamma_fibra_carbono
    area_rope = {'S': A_Small, 'M': A_Medium, 'L': A_Large}[Tamaño]
    mitad_barra = int(round((len(barras_grupo[0][0]) - 1) / 2))
    mitad_torre = int(round((len(torre) - 1) / 2))
    for esquina, barra_set in enumerate(barras_grupo):
        for barra in barra_set:
            nodo_mitad_barra = barra[mitad_barra][0]
            nodo_mitad_torre = torre[mitad_torre][esquina]
            ops.element('Truss', barra_actual, nodo_mitad_barra, nodo_mitad_torre, area_rope, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodo_mitad_barra, nodo_mitad_torre, Tamaño])
            barra_actual += 1
            nodo_final_barra = barra[-1][0]
            nodo_final_torre = torre[-1][esquina]
            ops.element('Truss', barra_actual, nodo_final_barra, nodo_final_torre, area_rope, 1, '-rho', gamma_fibra_carbono)
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
    caras_caja = [
        [1,4,8,5],
        [2,3,7,6],
        [1,2,6,5],
        [3,4,8,7],
        [1,2,3,4],
        [5,6,7,8]
    ]
    nodes_apoyos_list = [
        np.array([[3.3, 0, 1], [3.3, 2, -1], [3.3, -2, -1]]),
        np.array([[3.3, 3.9, 1], [0.3, 3.9, -1], [3.3, 0.9, -1]]),
        np.array([[0, 3.9, 1], [-2, 3.9, -1], [2, 3.9, -1]]),
        np.array([[-3.3, 3.9, 1], [-3.3, 0.9, -1], [-0.3, 3.9, -1]]),
        np.array([[-3.3, 0, 1], [-3.3, -2, -1], [-3.3, 2, -1]]),
        np.array([[-3.3, -3.9, 1], [-0.3, -3.9, -1], [-3.3, -0.9, -1]]),
        np.array([[0, -3.9, 1], [2, -3.9, -1], [-2, -3.9, -1]]),
        np.array([[3.3, -3.9, 1], [3.3, -0.9, -1], [0.3, -3.9, -1]])
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
    # Conexión de la torre con las barras
    barras_grupo_1 = [nodos_barras[7], nodos_barras[0], nodos_barras[1]]
    barras_grupo_2 = [nodos_barras[1], nodos_barras[2], nodos_barras[3]]
    barras_grupo_3 = [nodos_barras[3], nodos_barras[4], nodos_barras[5]]
    barras_grupo_4 = [nodos_barras[5], nodos_barras[6], nodos_barras[7]]
    conectar_mitades_torre_con_barras(nodos_torre_1, [barras_grupo_1, barras_grupo_2, barras_grupo_3, barras_grupo_4], Tamaño='L')


    # Acumulación de masas y análisis modal
    acumular_masa_barras_en_nodos()
    masa_barras = sum(masa_acumulada_por_nodo.values())
    area_paneles, masa_paneles = acumular_masa_paneles()
    print(f'Los paneles tienen un área total de {area_paneles} y una masa total de {masa_paneles}')
    print(f'La masa total de las barras es: {masa_barras}')
    print(f'La proporción es: {(masa_barras/masa_paneles)}')
    asignar_masas_a_nodos()
    eigenfrequencies = realizar_analisis_frecuencias(num_modes=10)
    visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

if __name__ == "__main__":
    main()
