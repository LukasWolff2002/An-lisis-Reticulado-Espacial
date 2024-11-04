import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import matplotlib.pyplot as plt
from variables import E, A, gamma, num_capas, unidad, gamma_rigido, masa_total, Area_panel, alpha, deltaT, barras_centrales

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
    return (nodo_i, nodo_j)

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

# --- Agregar Estructura Flexible ---
def agregar_estructura_flexible(nodos_conectados, A_flexible, gamma_flexible, element_id, elements, material_id_flexible):
    for i in range(len(nodos_conectados)):
        j = (i + 1) % len(nodos_conectados)
        elements.append(
            crear_elemento_truss(element_id, nodos_conectados[i], nodos_conectados[j], A_flexible, material_id_flexible, gamma_flexible)
        )
        element_id += 1

    return element_id

# --- Patrón de Inercia ---
def inercia ():
    paternx = 1
    paterny = 2
    paternz = 3

    ax = [0.1 * 9.81]  # Aceleración en x como lista
    ay = [0.1 * 9.81]  # Aceleración en y como lista
    az = [0.1 * 9.81]  # Aceleración en z como lista

    dirx = 1
    diry = 2
    dirz = 3

    ops.timeSeries('Constant', paternx, '-values', *ax)
    ops.timeSeries('Constant', paterny, '-values', *ay)
    ops.timeSeries('Constant', paternz, '-values', *az)

    ops.pattern('UniformExcitation', paternx, dirx, '-accel', paternx)
    ops.pattern('UniformExcitation', paterny, diry, '-accel', paterny)
    ops.pattern('UniformExcitation', paternz, dirz, '-accel', paternz)

# --- Asignación de Masa de Paneles ---
def acumular_masa_paneles(conexiones_paneles, masa_total, masa_acumulada_por_nodo):
    """
    Agrega la masa de los paneles a los nodos específicos que están conectados a los paneles.
    """
    masa_por_panel = masa_total / 4

    for nodos_conectados in conexiones_paneles:
        for nodo in nodos_conectados:
            masa_acumulada_por_nodo[nodo] = masa_acumulada_por_nodo.get(nodo, 0) + masa_por_panel

    return masa_acumulada_por_nodo

# Para asignar la masa acumulada a los nodos
def asignar_masas_a_nodos(masa_acumulada_por_nodo):
    for nodo, masa in masa_acumulada_por_nodo.items():
        print(f"Nodo {nodo} - Masa asignada: {masa}")
        ops.mass(nodo, masa, masa, masa)

# --- Acumulación de Masa de Barras en Nodos ---
def acumular_masa_barras_en_nodos(elements, A, gamma):
    """
    Calcula la masa de cada barra y la distribuye entre los nodos de sus extremos.
    Retorna un diccionario con la masa acumulada por cada nodo.
    """
    masa_acumulada_por_nodo = {}

    for element in elements:
        nodo_i, nodo_j = element
        
        # Obtener las coordenadas de los nodos
        coord_i = np.array(ops.nodeCoord(nodo_i))
        coord_j = np.array(ops.nodeCoord(nodo_j))
        
        # Calcular la longitud de la barra
        longitud = np.linalg.norm(coord_j - coord_i)
        
        # Calcular la masa de la barra
        masa_barra = longitud * A * gamma
        
        # Distribuir la mitad de la masa a cada nodo
        masa_acumulada_por_nodo[nodo_i] = masa_acumulada_por_nodo.get(nodo_i, 0) + masa_barra / 2
        masa_acumulada_por_nodo[nodo_j] = masa_acumulada_por_nodo.get(nodo_j, 0) + masa_barra / 2

    return masa_acumulada_por_nodo


# --- Visualización de la Estructura ---
def visualizar_estructura_rigida(nodos_paneles, plotter, color='green'):
    for nodos_conectados in nodos_paneles:
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.7)

# --- Graficar Desplazamientos Inerciales ---
def graficar_desplazamientos_inercia(elements, conexiones_paneles):
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
        truss.lines = np.hstack([[2, e[0] - 1, e[1] - 1] for e in elements])

        truss_displaced = truss.copy()
        truss_displaced.points[:, i] += displacements[i]
        plotter.add_mesh(truss, color='blue', label="Estructura Original")
        plotter.add_mesh(truss_displaced, color=colors[i], label=f"Estructura Desplazada en {labels[i]}")
        visualizar_estructura_rigida(conexiones_paneles, plotter, color='gold')
        plotter.add_text(f"Desplazamiento por Inercia en {labels[i]}", position='upper_left', font_size=10)
        plotter.show_axes()
        plotter.show()

    '''
    # Gráfico combinado de desplazamientos en X, Y y Z
    plotter_combined = plotters[3]
    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[0] - 1, e[1] - 1] for e in elements])

    truss_displaced_combined = truss.copy()
    truss_displaced_combined.points += np.column_stack(displacements)  # Aplicar desplazamientos en todas las direcciones

    plotter_combined.add_mesh(truss, color='blue', label="Estructura Original")
    plotter_combined.add_mesh(truss_displaced_combined, color=colors[3], label="Estructura Desplazada Combinada")
    visualizar_estructura_rigida(conexiones_paneles, plotter_combined, color='gold')
    plotter_combined.add_text("Desplazamiento por Inercia Combinado", position='upper_left', font_size=10)
    plotter_combined.show_axes()
    plotter_combined.show()
    '''

def calcular_masa_total_barras(elements, A, gamma):
    masa_total = 0
    for element in elements:
        nodo_i, nodo_j = element
        coord_i, coord_j = np.array(ops.nodeCoord(nodo_i)), np.array(ops.nodeCoord(nodo_j))
        longitud = np.linalg.norm(coord_j - coord_i)
        masa_barra = longitud * A * gamma
        masa_total += masa_barra
    return masa_total

def encontrar_barra(nodo_i, nodo_j, elements):
    """
    Encuentra el número de la barra que conecta dos nodos específicos.
    Parámetros:
        - nodo_i: ID del primer nodo.
        - nodo_j: ID del segundo nodo.
        - elements: Lista de elementos, donde cada elemento es una tupla (element_id, nodo_i, nodo_j).
    Retorna:
        - El número de la barra que conecta los nodos, o None si no existe.
    """
    for index, (ni, nj) in enumerate(elements):
        # Comprobar si el elemento conecta los nodos, independientemente del orden
        if (ni == nodo_i and nj == nodo_j) or (ni == nodo_j and nj == nodo_i):
            return index + 1  # Devolver el número de barra (basado en índice + 1)
    return None  # Devolver None si no se encuentra la barra

# --- Deformación Térmica ---
def variacion_termica(nodos_paneles, elements):
    termalStrain = alpha * deltaT
    ops.uniaxialMaterial('InitStrainMaterial', 2, 1, termalStrain)

    for listas in nodos_paneles:
        for i in range(len(listas)):
            nodo_i, nodo_j = listas[i], listas[(i + 1) % len(listas)]
            tag_barra = encontrar_barra(nodo_i, nodo_j, elements)
            if tag_barra is not None:
                ops.remove('element', tag_barra)
                ops.element('Truss', tag_barra, nodo_i, nodo_j, A, 2, '-rho', gamma)
                
                #print(f"Material `InitStrainMaterial` con deformación térmica aplicada a la barra {tag_barra}")

def graficar_desplazamientos_termicos(elements, nodos_paneles, escala=1.0):
    plotter = pv.Plotter()
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    deformaciones = [np.array([ops.nodeDisp(tag)[i] for tag in ops.getNodeTags()]) for i in range(3)]

    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[0] - 1, e[1] - 1] for e in elements])

    truss_desplazado = truss.copy()
    for i in range(3):
        truss_desplazado.points[:, i] += deformaciones[i] * escala

    plotter.add_mesh(truss, color='blue', label="Estructura Original")
    plotter.add_mesh(truss_desplazado, color='orange', label="Estructura Desplazada Térmicamente")
    visualizar_estructura_rigida(nodos_paneles, plotter, color='gold')
    plotter.add_text("Desplazamiento Térmico", position='upper_left', font_size=10)
    plotter.show_axes()
    plotter.show()

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

    
    E_rigido, material_id_rigido, A_rigido = 1e-1000, 999, 0.01
    ops.uniaxialMaterial('Elastic', material_id_rigido, E_rigido)

    conexiones_paneles = []
    for i in range(num_capas):
        conexiones_paneles.append([2 + 4 * i, 3 + 4 * i, 7 + 4 * i, 6 + 4 * i])
    for nodos_conectados in conexiones_paneles:
        element_id = agregar_estructura_flexible(nodos_conectados, A_rigido, gamma_rigido, element_id, elements, material_id_rigido)

    # Acumular masas de las barras en cada nodo
    masa_acumulada_por_nodo = acumular_masa_barras_en_nodos(elements, A, gamma)
    
    # Agregar la masa de los paneles a los nodos específicos conectados a los paneles
    masa_acumulada_por_nodo = acumular_masa_paneles(conexiones_paneles, masa_total, masa_acumulada_por_nodo)
    
    # Asignar las masas calculadas a los nodos
    asignar_masas_a_nodos(masa_acumulada_por_nodo)

    #ops.printModel('node')

    inercia()

    ops.system("BandGen")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    ops.analyze(100, 0.01)

    graficar_desplazamientos_inercia(elements, conexiones_paneles)

    variacion_termica(conexiones_paneles, elements)
    graficar_desplazamientos_termicos(elements, conexiones_paneles)



if __name__ == "__main__":
    main()
