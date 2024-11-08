# Importación de módulos necesarios
import numpy as np
import openseespy.opensees as ops
import pyvista as pv
import matplotlib.pyplot as plt
from variables import (
    E, A, gamma, num_capas, unidad, gamma_rigido, masa_total,
    Area_panel, alpha, deltaT, barras_centrales, cruz
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
        #elements.append(
            #crear_elemento_truss(element_id, nodos_conectados[i], nodos_conectados[j], A_flexible, material_id_flexible, gamma_flexible)
        #)
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
def graficar_desplazamientos_combinados_inercia(elements, conexiones_paneles, escala, titulo):
 
    plotter = pv.Plotter(window_size=[1920, 1080])
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
    plotter.screenshot(filename='INFORME/GRAFICOS/'+titulo+' '+cruz, transparent_background=False)


def graficar_desplazamientos_termicos(elements, nodos_paneles, escala):
    plotter = pv.Plotter(window_size=[1920, 1080])
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
    plotter.screenshot(filename='INFORME/GRAFICOS/Desplazamientos Termicos'+' '+cruz, transparent_background=False)


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
                continue
    
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

def graficar_barras_coloreadas(elements, fuerzas_barras, titulo, nodos_paneles):

    plotter = pv.Plotter(window_size=[1920, 1080])
    
    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    
    # Crear lista de líneas y colores
    lines = []
    colors = []
    for element in elements:
        element_id, nodo_i, nodo_j = element
        lines.append([2, nodo_i - 1, nodo_j - 1])
        
        fuerza = fuerzas_barras.get(element_id, 0)
        if fuerza > 0:
            color = 'blue'  # Tensión
        elif fuerza < 0:
            color = 'red'  # Compresión
        else:
            color = 'black'  # Sin fuerza
        colors.append(color)

    truss = pv.PolyData(nodes)
    truss.lines = np.hstack(lines)
    
    # Añadir cada barra con su color correspondiente
    for i, line in enumerate(truss.lines.reshape(-1, 3)):
        start_idx = line[1]
        end_idx = line[2]
        points = nodes[[start_idx, end_idx]]
        segment = pv.Line(points[0], points[1])
        plotter.add_mesh(segment, color=colors[i], line_width=2)
    
    # Agregar etiquetas de números de nodos (opcional)
    node_tags = ops.getNodeTags()
    labels = [str(tag) for tag in node_tags]
    plotter.add_point_labels(nodes, labels, point_size=5, font_size=12, text_color='black', name='labels')
    visualizar_estructura_rigida(nodos_paneles, plotter, color='gold')
    plotter.add_text(titulo, position='upper_left', font_size=10)
    plotter.add_legend([("Barras Tensionadas", "blue"), ("Barras Comprimidas", "red"),("Barras Sin Esfuerzo", "black")])
    plotter.show_axes()
    plotter.show()
    plotter.screenshot(filename='INFORME/GRAFICOS/'+titulo+' '+cruz, transparent_background=False)

def graficar_barras_coloreadas_por_esfuerzo_maximo(elements, nodes, fuerzas_barras, titulo):

    plotter = pv.Plotter(window_size=[1920, 1080])

    # Asegurarse de que 'nodes' es un array de numpy
    nodes = np.array(nodes)

    f_max = max(fuerzas_barras.values())
    
    #Ahora creo un diccionario con las fuerzas normalizadas
    fuerzas_normalizadas = {}

    for key in fuerzas_barras:
        fuerzas_normalizadas[key] = fuerzas_barras[key]/(f_max)

    # Crear lista de líneas y colores
    lines = []
    colors = []
    for element in elements:
        element_id, nodo_i, nodo_j = element
        lines.append([2, nodo_i - 1, nodo_j - 1])
        
        fuerza_normalizada = fuerzas_normalizadas.get(element_id, 0)
        
        colors.append(255*fuerza_normalizada)

    truss = pv.PolyData(nodes)
    truss.lines = np.hstack(lines)

    # Añadir cada barra con su color correspondiente
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
    plotter.add_text(titulo, position='upper_left', font_size=10)
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
    plotter.screenshot(filename='INFORME/GRAFICOS/'+titulo+' '+cruz, transparent_background=False)

    
    
    


def calculo_inercial (masa_acumulada_por_nodo, elements, conexiones_paneles, solucion, nodos_paneles):
    
    aceleracion_magnitud = 0.1 * 9.81  # Magnitud de la aceleración
    aceleraciones = [[aceleracion_magnitud,0,0, 'Aceleracion X', 300],[0,aceleracion_magnitud,0, 'Axeleracion Y', 300],[0,0,aceleracion_magnitud, 'Aceleracion Z', 1500]]

    fuerzas_maximas = {}

    for a in aceleraciones:
        fuerzas_barras = {}
        
        # Aplicar cargas inerciales
        aplicar_cargas_inerciales(masa_acumulada_por_nodo, a[0], a[1], a[2], solucion)

        # Configurar y ejecutar el análisis
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.algorithm('Linear')
        ops.analysis('Static')
        ops.analyze(1)

        graficar_desplazamientos_combinados_inercia(elements, conexiones_paneles, escala=a[4], titulo="Desplazamientos Inerciales "+a[3])
        obtener_reacciones_en_apoyos([1, 2, 3, 4])

        for element in elements:
            element_id, nodo_i, nodo_j = element

            # Obtener la fuerza interna en la barra
            fuerzas = ops.eleResponse(element_id, 'axialForce')
            fuerza_axial = fuerzas[0]  # Fuerza axial en el nodo i

            # Almacenar la fuerza máxima para la barra
            if element_id not in fuerzas_barras or abs(fuerza_axial) > abs(fuerzas_barras[element_id]):
                fuerzas_barras[element_id] = fuerza_axial

            if fuerza_axial not in fuerzas_maximas:
                fuerzas_maximas[element_id] = fuerza_axial**2

            else:
                fuerzas_maximas[element_id] += fuerza_axial**2

        # Limpiar cargas y análisis para la siguiente iteración
        ops.remove('loadPattern', solucion)
        ops.wipeAnalysis()

        solucion += 1
        
        # Graficar la estructura con las barras coloreadas al final
        graficar_barras_coloreadas(elements, fuerzas_barras, "Esfuerzos Internos Máximos en las Barras "+a[3], nodos_paneles)

    #Saco la raiz cuadrada
    for key in fuerzas_maximas:
        fuerzas_maximas[key] = (fuerzas_maximas[key])**0.5

    return solucion, fuerzas_maximas

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
    
    solucion, fuerzas_maximas = calculo_inercial(masa_acumulada_por_nodo, elements, conexiones_paneles, solucion, conexiones_paneles)

    
    print('\nLas fuerzas en las barras son\n',fuerzas_maximas)
    print('\nLas fuerzas maximas en las barras son\n',max(fuerzas_maximas.values()),'\n Y ocurre en la barra \n',max(fuerzas_maximas, key=fuerzas_maximas.get),'la cual une los nodos',elements[max(fuerzas_maximas, key=fuerzas_maximas.get)][1],'y',elements[max(fuerzas_maximas, key=fuerzas_maximas.get)][2])



    nodes = np.array([ops.nodeCoord(tag) for tag in ops.getNodeTags()])
    elements_barras = []
    for eleTag in ops.getEleTags():
        eleNodes = ops.eleNodes(eleTag)
        elements_barras.append((eleTag, eleNodes[0], eleNodes[1]))
    graficar_barras_coloreadas_por_esfuerzo_maximo(elements_barras, nodes, fuerzas_maximas, titulo="Esfuerzos Internos Máximos en las Barras " + cruz)

    #-----------------------------------------------------
    #Analisis Termico
    #-----------------------------------------------------


    # Aplicar deformación térmica y realizar análisis
    variacion_termica(conexiones_paneles, elements, solucion)
    
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)

    graficar_desplazamientos_termicos(elements, conexiones_paneles, escala=100)

    fuerzas_barras = {}
    for element in elements:
            element_id, nodo_i, nodo_j = element

            # Obtener la fuerza interna en la barra
            fuerzas = ops.eleResponse(element_id, 'axialForce')
            fuerza_axial = fuerzas[0]  # Fuerza axial en el nodo i

            # Almacenar la fuerza máxima para la barra
            if element_id not in fuerzas_barras or abs(fuerza_axial) > abs(fuerzas_barras[element_id]):
                fuerzas_barras[element_id] = fuerza_axial

    graficar_barras_coloreadas(elements, fuerzas_barras, "Esfuerzos Internos en las Barras por Variación Térmica", conexiones_paneles)

    masa_estructura = 0
    for masas in masa_acumulada_por_nodo.values():
        masa_estructura += masas

    print(f"\nMasa total de la estructura: {masa_estructura:.2f} kg")

if __name__ == "__main__":
    main()


#Debo encontrar theta y phi que maximicen la fuerza en cada barra
#Logicamente las barras que tienen mayor fuerz son las que unen el primer bloque

#la fuerza es maxima en cada barra cuando la acelaeracion es en la direccion de esta
