# Importación de módulos necesarios
import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import matplotlib.pyplot as plt
from variables import (
    E, A, gamma, num_capas, unidad, gamma_rigido, masa_total,
    Area_panel, alpha, deltaT, barras_centrales
)

# --- Configuración Inicial del Modelo ---
def inicializar_modelo():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

def definir_material(material_id, E):
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)

# --- Definición de Nodos y Elementos ---
def definir_nodos_iniciales(nodes):
    for i, n in enumerate(nodes):
        ops.node(i + 1, float(n[0]), float(n[1]), float(n[2]))

def definir_restricciones(nodo_tags, restricciones):
    for nodo in nodo_tags:
        ops.fix(nodo, *restricciones)

def crear_elemento_truss(element_id, nodo_i, nodo_j, A, material_id, gamma):
    if A <= 0:
        raise ValueError("El área de la sección transversal (A) debe ser mayor que 0.")
    if gamma < 0:
        raise ValueError("La densidad (gamma) no puede ser negativa.")
    ops.element('Truss', element_id, nodo_i, nodo_j, A, material_id, '-rho', gamma)
    return (element_id, nodo_i, nodo_j)

def conectar_nodos_cuadrado(node_tags, element_id, A, material_id, gamma, elements):
    for i in range(len(node_tags)):
        j = (i + 1) % len(node_tags)
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j], A, material_id, gamma))
        element_id += 1

    if node_tags[0] == 1:
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j]+1, A, material_id, gamma))
        element_id += 1

    elif node_tags[-1] == num_capas * 4 + 4:
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j]+1, A, material_id, gamma))
        element_id += 1

    else:
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j]+1, A, material_id, gamma))
        element_id += 1

        if barras_centrales:
            elements.append(crear_elemento_truss(element_id, node_tags[i]-1, node_tags[j], A, material_id, gamma))
            element_id += 1

    return element_id

def generar_capas(nodes, num_capas, element_id, A, material_id, gamma, elements):
    old_nodes = ops.getNodeTags()
    for capa in range(num_capas):
        last_node = ops.getNodeTags()[-1]
        new_nodes = []

        for j, n in enumerate(nodes):
            ops.node(j + last_node + 1, float(n[0]), float(n[1]), float(n[2] + unidad + (capa * unidad)))
            new_nodes.append(j + last_node + 1)

        element_id = conectar_nodos_cuadrado(new_nodes, element_id, A, material_id, gamma, elements)

        for i in old_nodes:
            j = new_nodes[0] if i == old_nodes[-1] else i + 5
            elements.append(crear_elemento_truss(element_id, i, j, A, material_id, gamma))
            element_id += 1

        for i in old_nodes:
            j = i + 4
            elements.append(crear_elemento_truss(element_id, i, j, A, material_id, gamma))
            element_id += 1

        old_nodes = new_nodes.copy()

    return element_id

# --- Agregar Estructura de Paneles ---
def agregar_estructura_flexible(nodos_conectados, A_flexible, gamma_flexible, element_id, elements, material_id_flexible):
    for i in range(len(nodos_conectados)):
        j = (i + 1) % len(nodos_conectados)
        elements.append(
            crear_elemento_truss(element_id, nodos_conectados[i], nodos_conectados[j], A_flexible, material_id_flexible, gamma_flexible)
        )
        element_id += 1
    return element_id

# --- Aplicación de Cargas Inerciales ---
def aplicar_cargas_inerciales(masa_acumulada_por_nodo):
    acceleration_x = 0.1 * 9.81  # Aceleración en X (m/s²)
    acceleration_y = 0.1 * 9.81  # Aceleración en Y (m/s²)
    acceleration_z = 0.1 * 9.81  # Aceleración en Z (m/s²)

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for nodo, masa in masa_acumulada_por_nodo.items():
        F_inertial_x = -masa * acceleration_x  # Fuerza inercial en X
        F_inertial_y = -masa * acceleration_y  # Fuerza inercial en Y
        F_inertial_z = -masa * acceleration_z  # Fuerza inercial en Z
        ops.load(nodo, F_inertial_x, F_inertial_y, F_inertial_z)

# --- Asignación de Masas ---
def acumular_masa_barras_en_nodos(elements, A, gamma):
    masa_acumulada_por_nodo = {}

    for element in elements:
        element_id, nodo_i, nodo_j = element
        coord_i = np.array(ops.nodeCoord(nodo_i))
        coord_j = np.array(ops.nodeCoord(nodo_j))
        longitud = np.linalg.norm(coord_j - coord_i)
        masa_barra = longitud * A * gamma
        masa_acumulada_por_nodo[nodo_i] = masa_acumulada_por_nodo.get(nodo_i, 0) + masa_barra / 2
        masa_acumulada_por_nodo[nodo_j] = masa_acumulada_por_nodo.get(nodo_j, 0) + masa_barra / 2

    return masa_acumulada_por_nodo

def acumular_masa_paneles(conexiones_paneles, masa_total, masa_acumulada_por_nodo):
    num_paneles = len(conexiones_paneles)
    masa_por_panel = masa_total / num_paneles

    nodos_paneles = []
    for nodos_conectados in conexiones_paneles:
        for nodo in nodos_conectados:
            masa_acumulada_por_nodo[nodo] = masa_acumulada_por_nodo.get(nodo, 0) + masa_por_panel / len(nodos_conectados)
            nodos_paneles.append(nodo)

    return masa_acumulada_por_nodo, nodos_paneles

def asignar_masas_a_nodos(masa_acumulada_por_nodo, nodos_paneles):
    for nodo, masa in masa_acumulada_por_nodo.items():
        if nodo not in nodos_paneles:
            print(f"Nodo {nodo} - Masa asignada (Barras): {masa}")
        else:
            print(f"Nodo {nodo} - Masa asignada (Barras + Panel): {masa}")
        ops.mass(nodo, masa, masa, masa)

# --- Visualización de la Estructura ---
def visualizar_estructura_rigida(nodos_paneles, plotter, color='green'):
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.7)

# --- Graficar Desplazamientos ---
def graficar_desplazamientos_inercia(elements, conexiones_paneles, escala):
    plotters = [pv.Plotter() for _ in range(3)]
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    displacements = [
        np.array([ops.nodeDisp(tag)[i] for tag in ops.getNodeTags()]) for i in range(3)
    ]

    print("Desplazamientos en X:", displacements[0])
    print("Desplazamientos en Y:", displacements[1])
    print("Desplazamientos en Z:", displacements[2])

    colors = ['cyan', 'green', 'red']
    labels = ['X', 'Y', 'Z']
    for i, plotter in enumerate(plotters):
        truss = pv.PolyData(nodes)
        truss.lines = np.hstack([[2, e[1] - 1, e[2] - 1] for e in elements])

        truss_displaced = truss.copy()
        truss_displaced.points[:, i] += displacements[i] * escala
        plotter.add_mesh(truss, color='blue', label="Estructura Original")
        plotter.add_mesh(truss_displaced, color=colors[i], label=f"Estructura Desplazada en {labels[i]}")
        visualizar_estructura_rigida(conexiones_paneles, plotter, color='gold')
        plotter.add_text(f"Desplazamiento por Inercia en {labels[i]} (Escala: {escala})", position='upper_left', font_size=10)
        plotter.show_axes()
        plotter.show()

def graficar_desplazamientos_termicos(elements, nodos_paneles, escala):
    plotter = pv.Plotter()
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    deformaciones = [np.array([ops.nodeDisp(tag)[i] for tag in ops.getNodeTags()]) for i in range(3)]

    print("Desplazamientos Termicos:", deformaciones[0])

    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[1] - 1, e[2] - 1] for e in elements])

    truss_desplazado = truss.copy()
    for i in range(3):
        truss_desplazado.points[:, i] += deformaciones[i] * escala

    plotter.add_mesh(truss, color='blue', label="Estructura Original")
    plotter.add_mesh(truss_desplazado, color='orange', label="Estructura Desplazada Térmicamente")
    visualizar_estructura_rigida(nodos_paneles, plotter, color='gold')
    plotter.add_text(f"Desplazamiento Térmico (Escala: {escala})", position='upper_left', font_size=10)
    plotter.show_axes()
    plotter.show()

# --- Funciones Auxiliares ---
def encontrar_barra(nodo_i, nodo_j, elements):
    for element in elements:
        element_id, ni, nj = element
        if (ni == nodo_i and nj == nodo_j) or (ni == nodo_j and nj == nodo_i):
            return element_id
    return None

# --- Aplicación de Deformación Térmica ---
def variacion_termica(conexiones_paneles, elements):
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    thermal_strain = alpha * deltaT

    for nodos_panel in conexiones_paneles:
        num_nodos = len(nodos_panel)
        for i in range(num_nodos):
            nodo_i = nodos_panel[i]
            nodo_j = nodos_panel[(i + 1) % num_nodos]
            element_id = encontrar_barra(nodo_i, nodo_j, elements)
            if element_id is not None:
                E_element = E
                A_element = A
                force_thermal = E_element * A_element * thermal_strain

                coord_i = np.array(ops.nodeCoord(nodo_i))
                coord_j = np.array(ops.nodeCoord(nodo_j))
                vector_unitario = (coord_j - coord_i) / np.linalg.norm(coord_j - coord_i)

                fuerza_nodo_i = -force_thermal * vector_unitario
                fuerza_nodo_j = force_thermal * vector_unitario

                ops.load(nodo_i, *fuerza_nodo_i)
                ops.load(nodo_j, *fuerza_nodo_j)
            else:
                print(f"No se encontró el elemento entre los nodos {nodo_i} y {nodo_j}")

# --- Ejecución del Análisis ---
def main():
    inicializar_modelo()
    definir_material(material_id=1, E=E)
    nodes = np.array([[0, 0, 0], [0, unidad, 0], [unidad, unidad, 0], [unidad, 0, 0]])
    definir_nodos_iniciales(nodes)
    definir_restricciones([1, 2, 3, 4], [1, 1, 1])

    elements = []
    element_id = 1
    element_id = conectar_nodos_cuadrado(ops.getNodeTags(), element_id, A, 1, gamma, elements)
    element_id = generar_capas(nodes, num_capas, element_id=element_id, A=A, material_id=1, gamma=gamma, elements=elements)

    # Definir material para los paneles
    E_Panel_Solar, material_id_rigido, A_rigido = 1e-1000, 999, 0.01
    ops.uniaxialMaterial('Elastic', material_id_rigido, E_Panel_Solar)

    # Definir conexiones de paneles
    conexiones_paneles = []
    for i in range(num_capas):
        conexiones_paneles.append([2 + 4 * i, 3 + 4 * i, 7 + 4 * i, 6 + 4 * i])
    for nodos_conectados in conexiones_paneles:
        element_id = agregar_estructura_flexible(nodos_conectados, A_rigido, gamma_rigido, element_id, elements, material_id_rigido)

    # Acumular y asignar masas
    masa_acumulada_por_nodo = acumular_masa_barras_en_nodos(elements, A, gamma)
    masa_acumulada_por_nodo, nodos_paneles = acumular_masa_paneles(conexiones_paneles, masa_total, masa_acumulada_por_nodo)
    asignar_masas_a_nodos(masa_acumulada_por_nodo, nodos_paneles)

    # Aplicar cargas inerciales y realizar análisis estático
    aplicar_cargas_inerciales(masa_acumulada_por_nodo)
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

    graficar_desplazamientos_inercia(elements, conexiones_paneles, escala=50)

    # Aplicar deformación térmica y realizar análisis
    variacion_termica(conexiones_paneles, elements)
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

    graficar_desplazamientos_termicos(elements, conexiones_paneles, escala=10)

if __name__ == "__main__":
    main()
