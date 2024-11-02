import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import matplotlib.pyplot as plt
from variables import E, A, gamma, num_capas, unidad, gamma_rigido, masa_total, Area_panel

#arch -x86_64 python3 /Users/lukaswolff/Desktop/24_20/METODOS_COMPUTACIONALES/Analisis-Reticulado-Espacial/CODIGO/ENTREGA_2/codigo.py

# Configuración inicial
def inicializar_modelo():
    ops.wipe()  # Limpiar cualquier modelo existente
    ops.model('basic', '-ndm', 3, '-ndf', 3)  # Modelo 3D con 3 grados de libertad

def definir_material(material_id, E):
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)

# Definición de nodos y elementos
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
        j = (i + 1) % len(node_tags)  # Cerrar el ciclo conectando el último con el primero
        #print(f'{i=}, {j=}')
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j], A, material_id, gamma))
        element_id += 1

    #Ahora conecto una diagonal
    if node_tags[0] == 1:
        #Estoy en el primer elemento
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j]+1, A, material_id, gamma))
        element_id += 1

    elif node_tags[-1] == num_capas*4 + 4:
        #Estoy en el ultimo elemento
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j]+1, A, material_id, gamma))
        element_id += 1

    else:
        #Agrego dos diagonales
        elements.append(crear_elemento_truss(element_id, node_tags[i], node_tags[j]+1, A, material_id, gamma))
        element_id += 1

        elements.append(crear_elemento_truss(element_id, node_tags[i]-1, node_tags[j], A, material_id, gamma))
        element_id += 1
    return element_id

# Generar múltiples capas de nodos y conectarlos
def generar_capas(nodes, num_capas, element_id, A, material_id, gamma, elements):
    old_nodes = ops.getNodeTags()  # Obtener los nodos existentes
    #print(num_capas)
    for capa in range(num_capas):
        #print(f'{capa=}')
        # Crear nueva capa de nodos desplazados
        last_node = ops.getNodeTags()[-1]
        new_nodes = []

        for j, n in enumerate(nodes):
            ops.node(j + last_node + 1, float(n[0]), float(n[1]), float(n[2] + unidad + (capa*unidad)))
            new_nodes.append(j + last_node + 1)

        # Conectar nodos nuevos entre sí
        element_id = conectar_nodos_cuadrado(new_nodes, element_id, A, material_id, gamma, elements)

        # Conectar diagonalmente entre capas
        for i in old_nodes:
            if i == old_nodes[-1]:
                j = new_nodes[0]
            else:
                j = i + 5
            elements.append(crear_elemento_truss(element_id, i, j, A, material_id, gamma))
            element_id += 1

        # Conectar verticalmente entre capas
        for i in old_nodes:
            j = i + 4
            elements.append(crear_elemento_truss(element_id, i, j, A, material_id, gamma))
            element_id += 1

        # Actualizar old_nodes para la siguiente capa
        old_nodes = new_nodes.copy()

    return element_id

# Agregar estructura rígida de 1m x 1m x 0.1m
def agregar_estructura_flexible(nodos_conectados, A_flexible, gamma_flexible, element_id, elements, material_id_flexible):
    """
    Agrega una estructura flexible conectada a ciertos nodos para evitar que influya en la rigidez global de la estructura.
    Parámetros:
        - nodos_conectados: Lista de nodos a los que se conectará la estructura flexible.
        - A_flexible: Área de la sección transversal de los elementos flexibles.
        - gamma_flexible: Densidad del material flexible.
        - element_id: ID inicial para los elementos.
        - elements: Lista donde se almacenan los elementos definidos.
        - material_id_flexible: Identificador del material flexible ya definido.
    """

    # Conectar los nodos con elementos flexibles
    for i in range(len(nodos_conectados)):
        j = (i + 1) % len(nodos_conectados)  # Cerrar el ciclo para crear el contorno
        elements.append(
            crear_elemento_truss(element_id, nodos_conectados[i], nodos_conectados[j], A_flexible, material_id_flexible, gamma_flexible)
        )
        element_id += 1

    return element_id



def acumular_masa_paneles(conexiones_paneles, masa_total):
    print(conexiones_paneles)
    """
    Acumula la masa total para cada nodo basado en la cantidad de paneles conectados.
    Parámetros:
        - conexiones_paneles: Lista de listas de nodos que forman cada panel.
        - masa_total: Masa total del panel (se distribuye uniformemente entre los nodos del panel).
    Retorna:
        - masa_acumulada_por_nodo: Diccionario con la masa acumulada para cada nodo.
    """
    masa_por_panel = masa_total / 4  # Cada panel tiene 4 nodos, la masa se distribuye uniformemente

    # Diccionario para almacenar la masa acumulada por nodo
    masa_acumulada_por_nodo = {}

    # Acumular la masa en cada nodo
    for nodos_conectados in conexiones_paneles:
        for nodo in nodos_conectados:
            if nodo not in masa_acumulada_por_nodo:
                masa_acumulada_por_nodo[nodo] = 0
            masa_acumulada_por_nodo[nodo] += masa_por_panel

    return masa_acumulada_por_nodo

def asignar_masas_a_nodos(masa_acumulada_por_nodo):
    """
    Asigna la masa acumulada a cada nodo en el modelo OpenSees.
    Parámetros:
        - masa_acumulada_por_nodo: Diccionario con la masa acumulada para cada nodo.
    """
    for nodo, masa in masa_acumulada_por_nodo.items():
        ops.mass(nodo, masa, masa, masa)  # Añadir masa en 3 grados de libertad


# Visualización con PyVista
def visualizar_estructura_rigida(nodos_paneles, plotter, color='green'):
    """
    Visualiza múltiples estructuras rígidas como paneles rellenos de color.
    Parámetros:
        - nodos_paneles: Lista de listas de nodos que forman cada panel.
        - plotter: Objeto PyVista para agregar la visualización.
        - color: Color de la superficie del panel.
    """
    for nodos_conectados in nodos_paneles:
        # Obtener las coordenadas de los nodos conectados
        puntos_panel = np.array([ops.nodeCoord(nodo) for nodo in nodos_conectados])

        # Crear una superficie plana (polygon) en PyVista
        surface = pv.PolyData(puntos_panel)
        surface.faces = [4, 0, 1, 2, 3]  # Definir una cara cuadrada con 4 vértices

        # Visualizar la superficie rellena de color
        plotter.add_mesh(surface, color=color, show_edges=True, opacity=0.7)

def visualizar_modo(mode_index, nodes, elements, eigenfrequencies, conexiones_paneles, scale=0.2):
    plotter = pv.Plotter()

    # Obtener la forma modal para el modo especificado
    #mode_shape = np.array([ops.nodeEigenvector(n, mode_index) for n in ops.getNodeTags()])
    
    #dl = np.sqrt(sum(((nodes[:, i].max() - nodes[:, i].min()) for i in range(3))))
    #max_mode = abs(mode_shape[:, :]).max()
    #scale_factor = scale * dl / max_mode

    # Aplicar desplazamientos escalados a las coordenadas del nodo
    #displaced_nodes = nodes + scale_factor * mode_shape

    # Convertir elements a formato 0-based para PyVista
    elements_0_based = np.array([(e[0] - 1, e[1] - 1) for e in elements])

    # Crear la malla del truss original
    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[0], e[1]] for e in elements_0_based])
    truss["mode"] = "original"

    # Crear la malla desplazada (deformada)
    #truss_deformed = pv.PolyData(displaced_nodes)
    #truss_deformed.lines = np.hstack([[2, e[0], e[1]] for e in elements_0_based])
    #truss_deformed["mode"] = f"Mode {mode_index}"

    # Configurar y mostrar el gráfico
    plotter.add_mesh(truss, color='blue', label="Original")
    #plotter.add_mesh(truss_deformed, color='red', label=f"Mode {mode_index}")
    
    # Agregar la visualización de todos los paneles rígidos
    visualizar_estructura_rigida(conexiones_paneles, plotter, color='gold')

    # Información de frecuencia
    #frequency = eigenfrequencies[mode_index - 1]
    #plotter.add_text(f"Mode {mode_index}: Frequency = {frequency:.4f} Hz", position='upper_left', font_size=10, color='black')
    #plotter.add_legend()
    #plotter.show()

# Graficar desplazamientos debido a la inercia en x, y, z, mostrando la estructura original y desplazada, incluyendo barras y paneles
def graficar_desplazamientos_inercia(elements, conexiones_paneles):
    # Crear plotters para cada dirección
    plotter_x = pv.Plotter()
    plotter_y = pv.Plotter()
    plotter_z = pv.Plotter()

    # Recuperar las coordenadas originales y desplazamientos nodales en cada dirección
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    displacements_x = np.array([ops.nodeDisp(tag)[0] for tag in ops.getNodeTags()])  # Desplazamiento en x
    displacements_y = np.array([ops.nodeDisp(tag)[1] for tag in ops.getNodeTags()])  # Desplazamiento en y
    displacements_z = np.array([ops.nodeDisp(tag)[2] for tag in ops.getNodeTags()])  # Desplazamiento en z

    # Crear la malla de la estructura original
    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[0] - 1, e[1] - 1] for e in elements])  # Crear conexiones para las barras

    # Agregar la estructura original y desplazada en cada dirección
    # Gráfico de desplazamiento en X
    truss_x = truss.copy()
    truss_x.points[:, 0] += displacements_x  # Aplicar desplazamiento en x
    plotter_x.add_mesh(truss, color='blue', label="Estructura Original")  # Estructura original
    plotter_x.add_mesh(truss_x, color='cyan', label="Estructura Desplazada en X")  # Estructura desplazada en X
    visualizar_estructura_rigida(conexiones_paneles, plotter_x, color='gold')  # Agregar paneles
    plotter_x.add_text("Desplazamiento por Inercia en X", position='upper_left', font_size=10)
    plotter_x.show()

    # Gráfico de desplazamiento en Y
    truss_y = truss.copy()
    truss_y.points[:, 1] += displacements_y  # Aplicar desplazamiento en y
    plotter_y.add_mesh(truss, color='blue', label="Estructura Original")  # Estructura original
    plotter_y.add_mesh(truss_y, color='green', label="Estructura Desplazada en Y")  # Estructura desplazada en Y
    visualizar_estructura_rigida(conexiones_paneles, plotter_y, color='gold')  # Agregar paneles
    plotter_y.add_text("Desplazamiento por Inercia en Y", position='upper_left', font_size=10)
    plotter_y.show()

    # Gráfico de desplazamiento en Z
    truss_z = truss.copy()
    truss_z.points[:, 2] += displacements_z  # Aplicar desplazamiento en z
    plotter_z.add_mesh(truss, color='blue', label="Estructura Original")  # Estructura original
    plotter_z.add_mesh(truss_z, color='red', label="Estructura Desplazada en Z")  # Estructura desplazada en Z
    visualizar_estructura_rigida(conexiones_paneles, plotter_z, color='gold')  # Agregar paneles
    plotter_z.add_text("Desplazamiento por Inercia en Z", position='upper_left', font_size=10)
    plotter_z.show()



# Realizar el análisis modal
def realizar_analisis_modal(num_modes):
    #Modicidar el tipo de analisis

    eigenvalues = ops.eigen(num_modes)
    eigenfrequencies = np.sqrt((eigenvalues)) / (2 * np.pi)
    print(f'Frecuencias naturales: {eigenfrequencies} Hz')
    return eigenfrequencies

def calcular_masa_total_barras(elements, A, gamma):
    """
    Calcula la masa total de todas las barras en el modelo.
    Parámetros:
        - elements: Lista de tuplas que contiene los pares de nodos que forman cada barra.
        - A: Área de la sección transversal de las barras.
        - gamma: Densidad del material de las barras.
    Retorna:
        - masa_total: Masa total de todas las barras en conjunto.
    """
    masa_total = 0

    for element in elements:
        nodo_i, nodo_j = element

        # Obtener coordenadas de los nodos
        coord_i = np.array(ops.nodeCoord(nodo_i))
        coord_j = np.array(ops.nodeCoord(nodo_j))

        # Calcular la longitud de la barra
        longitud = np.linalg.norm(coord_j - coord_i)

        # Calcular la masa de la barra
        masa_barra = longitud * A * gamma

        # Sumar la masa de la barra al total
        masa_total += masa_barra

    return masa_total

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



# Configuración y ejecución del análisis
def main():
    inicializar_modelo()
    definir_material(material_id=1, E=E)

    # Definir la geometría inicial del cuadrado
    nodes = np.array([[0, 0, 0], 
                      [0, unidad, 0], 
                      [unidad, unidad, 0], 
                      [unidad, 0, 0]])
    definir_nodos_iniciales(nodes)

    # Restricciones en los nodos iniciales
    definir_restricciones([1, 2, 3, 4], [1, 1, 1])

    # Crear los elementos iniciales y generar capas
    elements = []
    element_id = 1
    element_id = conectar_nodos_cuadrado(ops.getNodeTags(), element_id, A, 1, gamma, elements)
    element_id = generar_capas(nodes, num_capas, element_id=element_id, A=A, material_id=1, gamma=gamma, elements=elements)

    # Definir material rígido solo una vez
    E_rigido = 1e1  # Módulo de elasticidad alto para simular rigidez
    material_id_rigido = 999
    ops.uniaxialMaterial('Elastic', material_id_rigido, E_rigido)

    # Definir propiedades de la estructura rígida
    A_rigido = 0.01  # Área de sección transversal pequeña para simular el espesor

    # Nodos donde se conectará la estructura rígida
    nodos_totales = ops.getNodeTags()

    j = 2
    conexiones_paneles = [] 
    for i in range(num_capas):
        #print(j, j+1, j + 4, j+5)
        conexiones_paneles.append([j, j+1, j + 5, j+4])
        j += 4

    #Conexiones paneles es una lista con nodos a los cuales se les conecta un panel

    # Crear paneles
    for nodos_conectados in conexiones_paneles:
        element_id = agregar_estructura_flexible(nodos_conectados, A_rigido, gamma_rigido, element_id, elements, material_id_rigido)

    # Acumular masas en los nodos
    masa_acumulada_por_nodo = acumular_masa_paneles(conexiones_paneles, masa_total)
    
    # Asignar la masa acumulada a cada nodo
    asignar_masas_a_nodos(masa_acumulada_por_nodo)

    

    #visualizar_modo(i, np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()]), elements, 0, conexiones_paneles, scale=0.1)

    inercia()

    ops.system("BandGen")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    ops.analyze(100, 0.01)

    graficar_desplazamientos_inercia(elements, conexiones_paneles)



    """
    # Realizar análisis modal
    eigenfrequencies = realizar_analisis_modal(num_modes=10)

    # Visualizar modos de vibración
    for i in range(1, 11):
        visualizar_modo(i, np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()]), elements, eigenfrequencies, conexiones_paneles, scale=0.1)

    print('')
    masa_total_barras = calcular_masa_total_barras(elements, A, gamma)
    masa_total_panel = num_capas * gamma_rigido * Area_panel
    print(f'Masa total de todas las barras: {masa_total_barras} kg')

    print(f'Masa total de los paneles solares: {masa_total_panel} kg')

    masa_total_estructura = masa_total_barras + masa_total_panel

    porcentaje = (masa_total_barras/masa_total_panel)*100

    print(f'El porcentaje de masa es {porcentaje}%')

    ops.printModel("-node")
    """
    #ops.printModel("-node")
if __name__ == "__main__":
    main()
