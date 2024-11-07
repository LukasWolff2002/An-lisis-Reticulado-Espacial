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
def aplicar_cargas_inerciales(masa_acumulada_por_nodo, acceleration_x, acceleration_y, acceleration_z, solucion):

    ops.timeSeries('Constant', solucion)
    ops.pattern('Plain', solucion, solucion)
    for nodo, masa in masa_acumulada_por_nodo.items():
        F_inertial_x = masa * acceleration_x  # Fuerza inercial en X
        F_inertial_y = masa * acceleration_y  # Fuerza inercial en Y
        F_inertial_z = masa * acceleration_z  # Fuerza inercial en Z
        ops.load(nodo, F_inertial_x, F_inertial_y, F_inertial_z)

# --- Asignación de Masas ---
def acumular_masa_barras_en_nodos(elements, A, gamma):

    # Crear un diccionario para almacenar la masa acumulada por nodo
    masa_acumulada_por_nodo = {}

    # Obtener la lista de todos los nodos en el modelo
    todos_los_nodos = set()
    for element in elements:
        element_id, nodo_i, nodo_j = element
        todos_los_nodos.add(nodo_i)
        todos_los_nodos.add(nodo_j)

    # Para cada nodo, encontrar las barras conectadas a él
    for nodo in todos_los_nodos:
        masa_nodo = 0
        # Buscar elementos conectados al nodo actual
        for element in elements:
            element_id, nodo_i, nodo_j = element
            if nodo == nodo_i or nodo == nodo_j:
                # Calcular la masa de la barra
                coord_i = np.array(ops.nodeCoord(nodo_i))
                coord_j = np.array(ops.nodeCoord(nodo_j))
                longitud = np.linalg.norm(coord_j - coord_i)
                masa_barra = longitud * A * gamma
                # Agregar la mitad de la masa de la barra al nodo
                masa_nodo += masa_barra / 2
        # Almacenar la masa acumulada en el nodo
        masa_acumulada_por_nodo[nodo] = masa_nodo

    return masa_acumulada_por_nodo

def acumular_masa_paneles(conexiones_paneles, masa_total, masa_acumulada_por_nodo):
    num_paneles = len(conexiones_paneles)
    masa_por_panel = masa_total / 4

    nodos_paneles = []
    for nodos_conectados in conexiones_paneles:
        for nodo in nodos_conectados:
            #print(f"Nodo {nodo} - Masa asignada (Barras): {masa_acumulada_por_nodo.get(nodo, 0)}")
            masa_acumulada_por_nodo[nodo] = masa_acumulada_por_nodo.get(nodo, 0) + masa_por_panel
            nodos_paneles.append(nodo)

    return masa_acumulada_por_nodo, nodos_paneles

def asignar_masas_a_nodos(masa_acumulada_por_nodo, nodos_paneles):
    for nodo, masa in masa_acumulada_por_nodo.items():
        #if nodo not in nodos_paneles:
            #print(f"Nodo {nodo} - Masa asignada (Barras): {masa}")
       # else:
           # print(f"Nodo {nodo} - Masa asignada (Barras + Panel): {masa}")
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

    #print("Desplazamientos en X:", displacements[0])
    #print("Desplazamientos en Y:", displacements[1])
    #print("Desplazamientos en Z:", displacements[2])

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

def graficar_desplazamientos_combinados_inercia(elements, conexiones_paneles, escala=1.0, titulo="Desplazamiento Combinado"):
    """
    Grafica la estructura original y desplazada con desplazamientos combinados en X, Y y Z.
    
    Parámetros:
    - elements: Lista de elementos (element_id, nodo_i, nodo_j).
    - conexiones_paneles: Lista de paneles definidos por sus nodos.
    - escala: Factor de escala para amplificar los desplazamientos.
    - titulo: Título del gráfico.
    """
    plotter = pv.Plotter()
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    node_tags = ops.getNodeTags()
    displacements = [
        np.array([ops.nodeDisp(tag)[i] for tag in ops.getNodeTags()]) for i in range(3)
    ]
    #print("Desplazamientos en X:", displacements[0])
    #print("Desplazamientos en Y:", displacements[1])
    #print("Desplazamientos en Z:", displacements[2])

    displacements = np.array([ops.nodeDisp(tag) for tag in ops.getNodeTags()])

    #print("Desplazamientos Totales por Inercia:", displacements)
    
    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[1] - 1, e[2] - 1] for e in elements])
    
    truss_desplazado = truss.copy()
    truss_desplazado.points += displacements * escala
    
    plotter.add_mesh(truss, color='blue', label="Estructura Original")
    plotter.add_mesh(truss_desplazado, color='red', label="Estructura Desplazada")
    visualizar_estructura_rigida(conexiones_paneles, plotter, color='gold')
    plotter.add_legend([("Estructura Original", "blue"), ("Estructura Desplazada", "red")])
    plotter.add_text(f"{titulo} (Escala: {escala})", position='upper_left', font_size=10)
    
    # Agregar etiquetas de números de nodos
    labels = [str(tag) for tag in node_tags]
    plotter.add_point_labels(truss_desplazado.points, labels, point_size=5, font_size=12, text_color='black', name='labels')

    plotter.show_axes()
    plotter.show()


def graficar_desplazamientos_termicos(elements, nodos_paneles, escala):
    plotter = pv.Plotter()
    node_tags = ops.getNodeTags()
    nodes = np.array([ops.nodeCoord(tag) for tag in node_tags])
    deformaciones = [np.array([ops.nodeDisp(tag)[i] for tag in node_tags]) for i in range(3)]

    #print("Desplazamientos Térmicos X:", deformaciones[0])
    #print("Desplazamientos Térmicos Y:", deformaciones[1])
    #print("Desplazamientos Térmicos Z:", deformaciones[2])

    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[1] - 1, e[2] - 1] for e in elements])

    truss_desplazado = truss.copy()
    for i in range(3):
        truss_desplazado.points[:, i] += deformaciones[i] * escala

    plotter.add_mesh(truss, color='blue', label="Estructura Original")
    plotter.add_mesh(truss_desplazado, color='orange', label="Estructura Desplazada Térmicamente")
    visualizar_estructura_rigida(nodos_paneles, plotter, color='gold')
    plotter.add_text(f"Desplazamiento Térmico (Escala: {escala})", position='upper_left', font_size=10)

    # Agregar etiquetas de números de nodos
    labels = [str(tag) for tag in node_tags]
    plotter.add_point_labels(truss_desplazado.points, labels, point_size=5, font_size=12, text_color='black', name='labels')

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
def variacion_termica(conexiones_paneles, elements, solucion):
    ops.timeSeries('Linear', solucion)
    ops.pattern('Plain', solucion, solucion)
    thermal_strain = alpha * deltaT

    #Genero una lista con los nodos previamente cargados
    nodos_previamente_cargados = set()
    for nodos_panel in conexiones_paneles:
        num_nodos = len(nodos_panel)
        for i in range(num_nodos):
            nodo_i = nodos_panel[i]
            nodo_j = nodos_panel[(i + 1) % num_nodos]

            element_id = encontrar_barra(nodo_i, nodo_j, elements)

            par_nodos = frozenset({nodo_i, nodo_j})

            if par_nodos in nodos_previamente_cargados:
                print('Ya se cargo la barra', nodo_j, nodo_i)
    
            elif element_id is not None:
                E_element = E
                A_element = A
                force_thermal = E_element * A_element * thermal_strain

                coord_i = np.array(ops.nodeCoord(nodo_i))
                coord_j = np.array(ops.nodeCoord(nodo_j))
                vector_unitario = (coord_j - coord_i) / np.linalg.norm(coord_j - coord_i)

                fuerza_nodo_i = -force_thermal * vector_unitario
                fuerza_nodo_j = force_thermal * vector_unitario

                #print('se cargaron los nodos', nodo_i, nodo_j, 'con una fuerza termica' , fuerza_nodo_i, fuerza_nodo_j)

                ops.load(nodo_i, *fuerza_nodo_i)
                ops.load(nodo_j, *fuerza_nodo_j)

                nodos_previamente_cargados.add(par_nodos)
            #else:
                #print(f"No se encontró el elemento entre los nodos {nodo_i} y {nodo_j}")

def obtener_reacciones_en_apoyos(nodos_soporte):
    print("\nReacciones en los apoyos:")
    ops.reactions()
    for nodo in nodos_soporte:
        reacciones = ops.nodeReaction(nodo)
        print(f"Nodo {nodo}: Reacción Rx = {reacciones[0]:.2f}, Ry = {reacciones[1]:.2f}, Rz = {reacciones[2]:.2f}")

def optimizador_inercial (masa_acumulada_por_nodo, elements, conexiones_paneles, solucion):
    # Las aceleraciones que se estan aplicando al analisis son:
    acceleration_x = 0.1 * 9.81  # Aceleración en X (m/s²)
    acceleration_y = 0.1 * 9.81  # Aceleración en Y (m/s²)
    acceleration_z = 0.1 * 9.81  # Aceleración en Z (m/s²)

    fuerzas_barras = {}

    aceleraciones = [[acceleration_x,0,0],[0, acceleration_y, 0],[0,0, acceleration_z]]

    for a in aceleraciones:

        aplicar_cargas_inerciales(masa_acumulada_por_nodo, a[0], a[1], a[2], solucion)
        
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.algorithm('Linear')
        ops.analysis('Static')
        ops.analyze(1)

        #Esta funcion grafica los desplzamientos en cada direccion por separado
        #graficar_desplazamientos_inercia(elements, conexiones_paneles, escala=50)

        #Esta funcion grafica los desplazamientos en las tres direcciones combinados
        graficar_desplazamientos_combinados_inercia(elements, conexiones_paneles, escala=100, titulo="Desplazamiento Inercial")

        # Obtener reacciones en los apoyos
        nodos_soporte = [1, 2, 3, 4]  # Asegúrate de que esta lista contenga los nodos correctos
        obtener_reacciones_en_apoyos(nodos_soporte)

        print("Fuerzas internas en cada barra:")
        
        for element in elements:
            element_id, nodo_i, nodo_j = element
            # Obtener la fuerza interna en el elemento
            # Para elementos Truss, 'force' devuelve [N_i, N_j], donde N_i = -N_j
            fuerzas = ops.eleResponse(element_id, 'axialForce')
            # La fuerza axial es igual en magnitud y opuesta en dirección en ambos nodos
            fuerza_axial = fuerzas[0]  # Tomamos la fuerza en el nodo i
            if element_id not in fuerzas_barras:
                fuerzas_barras[element_id] = fuerza_axial**2
            else:
                fuerzas_barras[element_id] += fuerza_axial**2
            
        
        # Resetear el estado del modelo antes del análisis térmico o para el siguiente analisis inercial
        #ops.remove('loadPattern', 1)  # Si el patrón de carga inercial tiene tag 1
        #ops.remove('loadPattern', 2)  # Si el patrón de carga inercial tiene tag 2
        #ops.remove('loadPattern', 3)  # Si el patrón de carga inercial tiene tag 3
        # Limpiar análisis anterior
        ops.remove('loadPattern', solucion)
        ops.wipeAnalysis()

        solucion += 1
    
    for barras in fuerzas_barras:
        fuerzas_barras[barras] = np.sqrt(fuerzas_barras[barras])

    return solucion, fuerzas_barras

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

    solucion = 1 #es el numero de analisi realizado para timeseries
    #-----------------------------------------------------
    #Analisis Inercial
    #-----------------------------------------------------
    
    solucion, fuerzas_barras = optimizador_inercial(masa_acumulada_por_nodo, elements, conexiones_paneles, solucion)

    
    print('\nLas fuerzas en las barras son\n',fuerzas_barras)

    #-----------------------------------------------------
    #Analisis Termico
    #-----------------------------------------------------


    # Aplicar deformación térmica y realizar análisis
    #variacion_termica(conexiones_paneles, elements, solucion)
    
    #ops.system('BandSPD')
    #ops.numberer('RCM')
    #ops.constraints('Plain')
    #ops.integrator('LoadControl', 1.0)
    #ops.algorithm('Linear')
    #ops.analysis('Static')
    #ops.analyze(1)

    #graficar_desplazamientos_termicos(elements, conexiones_paneles, escala=100)

    masa_estructura = 0
    for masas in masa_acumulada_por_nodo.values():
        masa_estructura += masas

    print(f"\nMasa total de la estructura: {masa_estructura:.2f} kg")

if __name__ == "__main__":
    main()


#Debo encontrar theta y phi que maximicen la fuerza en cada barra
#Logicamente las barras que tienen mayor fuerz son las que unen el primer bloque

#la fuerza es maxima en cada barra cuando la acelaeracion es en la direccion de esta
