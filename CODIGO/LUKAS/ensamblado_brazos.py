import numpy as np
import openseespy.opensees as ops
import pyvista as pv
from math import atan

barras = []
nodos_barra_1, nodos_barra_2, nodos_barra_3, nodos_barra_4, nodos_barra_5, nodos_barra_6, nodos_barra_7, nodos_barra_8 = [], [], [], [], [], [], [], []
largo_barras, espaciamiento = 10, 5

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

def nodos_apoyos (nodos_apoyo, nodos_barra, nodo_actual):
    for i, n in enumerate(nodos_apoyo):
        nodo_actual += 1
        ops.node(nodo_actual, float(n[0]), float(n[1]), float(n[2]))
        #Dejo fijo el nodo
        ops.fix(nodo_actual, 1, 1, 1)

    nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

def generar_elemento_axial(nodo_1, nodo_2):
    global barra_actual, A_Small, gamma_fibra_carbono
    ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Small, 1, '-rho', gamma_fibra_carbono)
    barras.append([barra_actual, nodo_1, nodo_2])
    barra_actual += 1

def nodos_base_barras(nodos_iniciales_barra, nodos_barra, alpha):
    global nodo_actual, barra_actual, ancho_barras, alto_barras, largo_inicial_barras
    nodos_iniciales_barra = [largo_inicial_barras, alpha, 0]
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
        generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j])

    #Ahora conecto con los apoyos
    if len(nodos_barra) > 1:
        conectar_capas(nodos_barra[-2], nodos_barra[-1])

def conectar_capas (nodos_capa_1, nodos_capa_2):
    global barra_actual, A_Small, gamma_fibra_carbono
    #Hago la conexión entre las capas
    for i in range(len(nodos_capa_1)):
        j = (i + 1) % len(nodos_capa_1)
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[i], A_Small, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[i]])
        barra_actual += 1

    #Agrego la conexion diagonal
    for i in range(len(nodos_capa_1)):
        ops.element('Truss', barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)], A_Small, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodos_capa_1[i], nodos_capa_2[(i+1)%len(nodos_capa_2)]])
        barra_actual += 1

def alargar_barras (alpha, nodos_barra):
    global nodo_actual
    for i in range(int(largo_barras/espaciamiento)):
        largo_adicional = (i+1)*espaciamiento
        nodos_nuevos = [largo_inicial_barras + largo_adicional, alpha, 0]
        x, y, z = coordenadas_cartesianas(nodos_nuevos[0], nodos_nuevos[1], nodos_nuevos[2])
        x_1, y_1 = calcular_nuevo_punto_transversal(nodos_nuevos[1], nodos_nuevos[0], ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(nodos_nuevos[1], nodos_nuevos[0], -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z - alto_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z - alto_barras)

        nodos_barra.append([nodo_actual - 2, nodo_actual - 1, nodo_actual])

        for i in range(len(nodos_barra[-1])):
            j = (i + 1) % len(nodos_barra[-1])
            generar_elemento_axial(nodos_barra[-1][i], nodos_barra[-1][j])

        conectar_capas(nodos_barra[-2], nodos_barra[-1])


def ensamblar_brazos(nodo_actual):
    
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
    
    nodos_apoyos(nodes_apoy_1, nodos_barra_1, nodo_actual)
    nodos_apoyos(nodes_apoy_2, nodos_barra_2, nodo_actual)
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
