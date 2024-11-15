import numpy as np
import openseespy.opensees as ops
import pyvista as pv
from math import atan

nodo_actual, barra_actual = 0, 0

def inicializar_modelo():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

def definir_material(material_id, E):
    if E <= 0:
        raise ValueError("El módulo de elasticidad debe ser mayor que 0.")
    ops.uniaxialMaterial('Elastic', material_id, E)

def definir_nodos(nodes):
    for i, n in enumerate(nodes):
        ops.node(i + 1, float(n[0]), float(n[1]), float(n[2]))

def definir_seccion_tubo(seccion_id, D1, D2):
    if D1 <= D2:
        raise ValueError("El diámetro exterior debe ser mayor que el diámetro interior.")
    A = np.pi * (D1**2 - D2**2) / 4
    return A

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
    plotter.add_point_labels(points, labels, font_size=30, text_color="red", always_visible=True)

    # Visualizar las caras de la caja satélite como paneles
    visualizar_panel(caras_caja, plotter, color='blue')
    visualizar_panel(conexiones_paneles, plotter, color='gold')

    # Dibujar las barras que conectan los nodos
    for barra in barras:
        _, nodo_i, nodo_j = barra  # Ignorar el número de barra, solo usar nodo_i y nodo_j
        nodo_inicio, nodo_fin = nodo_i - 1, nodo_j - 1  # Restar 1 para índices de numpy
        coord_inicio = nodes[nodo_inicio]
        coord_fin = nodes[nodo_fin]
        linea = pv.Line(coord_inicio, coord_fin)
        plotter.add_mesh(linea, color="black", line_width=2)  # Agregar barra entre los nodos

    # Mostrar ejes y renderizar la visualización
    plotter.show_axes()
    plotter.show()



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

def angulos_barras ():

        alpha_1 = np.degrees(atan(1.95/1.65))
        alpha_2 = alpha_1 + np.degrees(2*atan(1.65/1.95))
        alpha_3 = 2*alpha_1 + alpha_2
        alpha_4 = alpha_3 + np.degrees(2*atan(1.65/1.95))

        #Aqui puedo mover las barras
        alpha_1 += 0
        alpha_2 += 0
        alpha_3 += 0
        alpha_4 += 0

        return alpha_1, alpha_2, alpha_3, alpha_4

def generar_elemento_axial (nodo_1, nodo_2):
        global nodo_actual, barra_actual, A_Small, gamma_fibra_carbono
        ops.element('Truss', barra_actual, nodo_1, nodo_2, A_Small, 1, '-rho', gamma_fibra_carbono)
        barras.append([barra_actual, nodo_1, nodo_2])
        barra_actual += 1

#Ahora genero la primera seccion se nodos, de las cuales se alargara cada barra
def nodos_base_barras (nodos_iniciales_barra, nodos_barra, alpha):
    global nodo_actual, barra_actual
    nodos_iniciales_barra = [largo_inicial_barras, alpha, 0]

    x, y, z = coordenadas_cartesianas(nodos_iniciales_barra[0], nodos_iniciales_barra[1], nodos_iniciales_barra[2])
    x_1, y_1 = calcular_nuevo_punto_transversal(nodos_iniciales_barra[1], nodos_iniciales_barra[0], ancho_barras)
    x_2, y_2 = calcular_nuevo_punto_transversal(nodos_iniciales_barra[1], nodos_iniciales_barra[0], -ancho_barras)

    nodo_actual += 1
    ops.node(nodo_actual, x, y, z)
    nodo_actual += 1
    ops.node(nodo_actual, x_1, y_1, z-alto_barras)
    nodo_actual += 1
    ops.node(nodo_actual, x_2, y_2, z-alto_barras)

    nodos_barra.append([nodo_actual-2, nodo_actual-1, nodo_actual])

    #Uno los 3 nodos
    for i in range(len(nodos_barra[0])):
        j = ((i + 1) % len(range(3))) 
        generar_elemento_axial(nodos_barra[0][i], nodos_barra[0][j])

def main():
    global nodo_actual, barra_actual
    # -----------------------------------
    # Definición del modelo
    # -----------------------------------
    
    # Se inicia un modelo con 3 grados de libertad por nodo
    inicializar_modelo()

    # -----------------------------------
    # Definición de materiales y tipos de barras
    # -----------------------------------

    # Defino el material fibra de carbono
    gamma_fibra_carbono = 1.91*1000
    E_fibra_carbono = 338e9 
    definir_material(material_id=1, E=E_fibra_carbono)

    # Defino barra pequeña
    D1_Small = 0.0058
    D2_Small = 0.002400
    A_Small = definir_seccion_tubo(seccion_id=1, D1=D1_Small, D2=D2_Small)

    # Defino barra mediana
    D1_Medium = 0.008
    D2_Medium = 0.004
    A_Medium = definir_seccion_tubo(seccion_id=2, D1=D1_Medium, D2=D2_Medium)

    # Defino barra grande
    D1_Large = 0.010
    D2_Large = 0.006
    A_Large = definir_seccion_tubo(seccion_id=3, D1=D1_Large, D2=D2_Large)

    # -----------------------------------
    # Genero la caja satélite
    # -----------------------------------

    # Primero compongo los 8 nodos
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

    nodo_actual += len(nodes)


    # Defino las caras, ajustando los índices para que coincidan con los índices de OpenSees (iniciando desde 1)
    caras_caja = [
        [1,4,8,5],
        [2,3,7,6],
        [1,2,6,5],
        [3,4,8,7],
        [1,2,3,4],
        [5,6,7,8]
    ]

    # Genero los nodos
    definir_nodos(nodes)

    #-----------------------------------
    #Genero primero elemento de cada barra
    #-----------------------------------

    #Defino el angulo que va a tener cada barra
    global alpha_1, alpha_2, alpha_3, alpha_4
    alpha_1, alpha_2, alpha_3, alpha_4 = angulos_barras()

    
    ancho_barras = 1 #metro
    alto_barras = 1 #metros
    largo_inicial_barras = 7 #metros

    nodos_barra_1 = []
    nodos_barra_2 = []
    nodos_barra_3 = []
    nodos_barra_4 = []

    barras = []

    #[largo, angulo_xy, altura_z]
    nodes_iniciales_barras_1 = [largo_inicial_barras, alpha_1, 0]
    nodes_iniciales_barras_2 = [largo_inicial_barras, alpha_2, 0]
    nodes_iniciales_barras_3 = [largo_inicial_barras, alpha_3, 0]
    nodes_iniciales_barras_4 = [largo_inicial_barras, alpha_4, 0]

    print(nodo_actual)

    nodos_base_barras(nodes_iniciales_barras_1, nodos_barra_1, alpha_1)
    nodos_base_barras(nodes_iniciales_barras_2, nodos_barra_2, alpha_2)
    nodos_base_barras(nodes_iniciales_barras_3, nodos_barra_3, alpha_3)
    nodos_base_barras(nodes_iniciales_barras_4, nodos_barra_4, alpha_4)


    #El formato es

    
    

    

    largo_barras = 10 #metros
    espaciamiento = 5 #metros

    for i in range(int(largo_barras/espaciamiento)):

        largo_adicional = (i+1)*espaciamiento
        nodes_nuevos_barras_1 = [largo_inicial_barras+largo_adicional, alpha_1, 0]
        x, y, z = coordenadas_cartesianas(nodes_nuevos_barras_1[0], nodes_nuevos_barras_1[1], nodes_nuevos_barras_1[2])
        x_1, y_1 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_1[1], nodes_nuevos_barras_1[0], ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_1[1], nodes_nuevos_barras_1[0], -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z-alto_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z-alto_barras)

        nodos_barra_1.append([nodo_actual-2, nodo_actual-1, nodo_actual])

        #Ahora uno los 3 nodos

        for i in range(len(nodos_barra_1[-1])):
            j = ((i + 1) % len(range(3))) 

            print(nodos_barra_1[-1][i], nodos_barra_1[-1][j])
            
            #por lo tanto, unos los nodos
            ops.element('Truss', barra_actual, nodos_barra_1[-1][i], nodos_barra_1[-1][j], A_Small, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_1[-1][i], nodos_barra_1[-1][j]])
            barra_actual += 1

        #Ahora uno las capas
        for i in range(len(nodos_barra_1[-2])):
            
            ops.element('Truss', barra_actual, nodos_barra_1[-2][i], nodos_barra_1[-1][i], A_Small, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_1[-2][i], nodos_barra_1[-1][i]])
            barra_actual += 1

        #AHora uno las capas con diagonales
        for i in range(len(nodos_barra_1[-2])):
            j = ((i + 1) % len(range(3))) 

            ops.element('Truss', barra_actual, nodos_barra_1[-2][i], nodos_barra_1[-1][j], A_Small, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_1[-2][i], nodos_barra_1[-1][j]])
            barra_actual += 1


        nodes_nuevos_barras_2 = [largo_inicial_barras+largo_adicional, alpha_2, 0]
        x, y, z = coordenadas_cartesianas(nodes_nuevos_barras_2[0], nodes_nuevos_barras_2[1], nodes_nuevos_barras_2[2])
        x_1, y_1 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_2[1], nodes_nuevos_barras_2[0], ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_2[1], nodes_nuevos_barras_2[0], -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z-alto_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z-alto_barras)

        nodos_barra_2.append([nodo_actual-2, nodo_actual-1, nodo_actual])

        #Ahora uno los 3 nodos
        for i in range(len(nodos_barra_2[-1])):
            j = ((i + 1) % len(range(3))) 

            print(nodos_barra_2[-1][i], nodos_barra_2[-1][j])
            
            #por lo tanto, unos los nodos
            ops.element('Truss', barra_actual, nodos_barra_2[-1][i], nodos_barra_2[-1][j], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_2[-1][i], nodos_barra_2[-1][j]])
            barra_actual += 1

        #Ahora uno las capas
        for i in range(len(nodos_barra_2[-2])):
                
            ops.element('Truss', barra_actual, nodos_barra_2[-2][i], nodos_barra_2[-1][i], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_2[-2][i], nodos_barra_2[-1][i]])
            barra_actual += 1

        #AHora uno las capas con diagonales
        for i in range(len(nodos_barra_2[-2])):
            j = ((i + 1) % len(range(3))) 

            ops.element('Truss', barra_actual, nodos_barra_2[-2][i], nodos_barra_2[-1][j], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_2[-2][i], nodos_barra_2[-1][j]])
            barra_actual += 1

        nodes_nuevos_barras_3 = [largo_inicial_barras+largo_adicional, alpha_3, 0]
        x, y, z = coordenadas_cartesianas(nodes_nuevos_barras_3[0], nodes_nuevos_barras_3[1], nodes_nuevos_barras_3[2])
        x_1, y_1 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_3[1], nodes_nuevos_barras_3[0], ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_3[1], nodes_nuevos_barras_3[0], -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z-alto_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z-alto_barras)

        nodos_barra_3.append([nodo_actual-2, nodo_actual-1, nodo_actual])

        #Ahora uno los 3 nodos

        for i in range(len(nodos_barra_3[-1])):
            j = ((i + 1) % len(range(3))) 

            print(nodos_barra_3[-1][i], nodos_barra_3[-1][j])
            
            #por lo tanto, unos los nodos
            ops.element('Truss', barra_actual, nodos_barra_3[-1][i], nodos_barra_3[-1][j], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_3[-1][i], nodos_barra_3[-1][j]])
            barra_actual += 1

        #Ahora uno las capas
        for i in range(len(nodos_barra_3[-2])):
                
            ops.element('Truss', barra_actual, nodos_barra_3[-2][i], nodos_barra_3[-1][i], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_3[-2][i], nodos_barra_3[-1][i]])
            barra_actual += 1

        #AHora uno las capas con diagonales

        for i in range(len(nodos_barra_3[-2])):
            j = ((i + 1) % len(range(3))) 

            ops.element('Truss', barra_actual, nodos_barra_3[-2][i], nodos_barra_3[-1][j], A_Medium, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_3[-2][i], nodos_barra_3[-1][j]])
            barra_actual += 1

        nodes_nuevos_barras_4 = [largo_inicial_barras+largo_adicional, alpha_4, 0]
        x, y, z = coordenadas_cartesianas(nodes_nuevos_barras_4[0], nodes_nuevos_barras_4[1], nodes_nuevos_barras_4[2])
        x_1, y_1 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_4[1], nodes_nuevos_barras_4[0], ancho_barras)
        x_2, y_2 = calcular_nuevo_punto_transversal(nodes_nuevos_barras_4[1], nodes_nuevos_barras_4[0], -ancho_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x, y, z)
        nodo_actual += 1
        ops.node(nodo_actual, x_1, y_1, z-alto_barras)
        nodo_actual += 1
        ops.node(nodo_actual, x_2, y_2, z-alto_barras)

        nodos_barra_4.append([nodo_actual-2, nodo_actual-1, nodo_actual])

        #Ahora uno los 3 nodos

        for i in range(len(nodos_barra_4[-1])):

            j = ((i + 1) % len(range(3))) 

            print(nodos_barra_4[-1][i], nodos_barra_4[-1][j])
            
            #por lo tanto, unos los nodos
            ops.element('Truss', barra_actual, nodos_barra_4[-1][i], nodos_barra_4[-1][j], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_4[-1][i], nodos_barra_4[-1][j]])
            barra_actual += 1

        #Ahora uno las capas

        for i in range(len(nodos_barra_4[-2])):

            ops.element('Truss', barra_actual, nodos_barra_4[-2][i], nodos_barra_4[-1][i], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_4[-2][i], nodos_barra_4[-1][i]])
            barra_actual += 1

        #AHora uno las capas con diagonales

        for i in range(len(nodos_barra_4[-2])):

            j = ((i + 1) % len(range(3))) 

            ops.element('Truss', barra_actual, nodos_barra_4[-2][i], nodos_barra_4[-1][j], A_Large, 1, '-rho', gamma_fibra_carbono)
            barras.append([barra_actual, nodos_barra_4[-2][i], nodos_barra_4[-1][j]])
            barra_actual += 1


    print(nodos_barra_1, nodos_barra_2, nodos_barra_3, nodos_barra_4)
    conexiones_paneles = []
    #Ahora creo los paneles solares
    for i in range((len(nodos_barra_1)-1)):
        
        #Entre las barras 1 y 2
        #print(nodos_barra_1[i][0], nodos_barra_2[i][0],nodos_barra_2[i+1][0],nodos_barra_1[i+1][0])
        conexiones_paneles.append([nodos_barra_1[i][0], nodos_barra_2[i][0],nodos_barra_2[i+1][0],nodos_barra_1[i+1][0]])

        #entre las barras 2 y 3
        #print(nodos_barra_2[i][0], nodos_barra_3[i][0],nodos_barra_3[i+1][0],nodos_barra_2[i+1][0])
        conexiones_paneles.append([nodos_barra_2[i][0], nodos_barra_3[i][0],nodos_barra_3[i+1][0],nodos_barra_2[i+1][0]])

        #entre las barras 3 y 4
        #print(nodos_barra_3[i][0], nodos_barra_4[i][0],nodos_barra_4[i+1][0],nodos_barra_3[i+1][0])
        conexiones_paneles.append([nodos_barra_3[i][0], nodos_barra_4[i][0],nodos_barra_4[i+1][0],nodos_barra_3[i+1][0]])

        #entre las barras 4 y 1
        #print(nodos_barra_4[i][0], nodos_barra_1[i][0],nodos_barra_1[i+1][0],nodos_barra_4[i+1][0])
        conexiones_paneles.append([nodos_barra_4[i][0], nodos_barra_1[i][0],nodos_barra_1[i+1][0],nodos_barra_4[i+1][0]])
    
    # Visualizo la caja satélite
    visualizar_caja_satelite(caras_caja, barras, conexiones_paneles)

if __name__ == "__main__":
    main()