import numpy as np
import math
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import pyvista as pv

from variables import E, A, gamma, num_capas, unidad

#arch -x86_64 python3 /Users/lukaswolff/Desktop/24_20/METODOS_COMPUTACIONALES/Analisis-Reticulado-Espacial/CODIGO/prueba_lukas_inicial.py

#Creo las barras con una gamma rho y se aplica directo al modelo

#La matriz de masa es una matriz diagonal que tiene la masa de cada elemento en la diagonal
#[m1, 0, 0]
#[0, m2, 0]
#[0, 0, m3]

#La masa es la del nodo

#La masa que del nodo es la mitad de todas las barras que le llegan
#Tambien  hay que considerar la masa de los paneles solares que se conectan a los nodos

#Van a haber nodos conectados al servidor, hay que ver que  porcion de masa toma cada uno

#Hay que definir la matriz de rigidez

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)


#opensees . eigen
#Aqui no vemos el desplzamiento por una fuerza, sino el desplazamiento por frecuencia

#Partamos definiendo un cuadrado y lo extendemos en 3D simetricamente

#Al agregar las masas a los modos, hay que definirlas en 3 grados de libertad


nodes = np.array([[0, 0, 0],
                  [0, unidad, 0],
                  [unidad, unidad, 0],
                  [unidad, 0, 0]])


# Definir nodos en OpenSees
for i, n in enumerate(nodes):
    #print(i,n)
    ops.node(i+1, float(n[0]), float(n[1]), float(n[2]))

#Defino los nodos de apoyo
ops.fix(1, 1, 1, 0)
ops.fix(2, 1, 1, 0)
ops.fix(3, 1, 1, 0)
ops.fix(4, 1, 1, 0)


# Definir el material uniaxial elástico
ops.uniaxialMaterial('Elastic', 1, E)

#Obtengo todos los nodos
node_tags = ops.getNodeTags()

# Crear una lista de elementos conectando los nodos
elements = []
element_id = 1  # Iniciar el identificador de elementos

for i in range(len(node_tags)):
    j = (i + 1) % len(node_tags)  # Usar el módulo para cerrar el ciclo
    elements.append((node_tags[i], node_tags[j]))  # Guardar las conexiones
    
    ops.element('Truss', element_id, node_tags[i], node_tags[j], A, 1, '-rho', gamma)
    element_id += 1  # Incrementar el identificador del elemento para que sea único

    #El element id no lo debo resetear nunca















old_nodes = ops.getNodeTags()
#Bien, ya con el cuadrado definido, vamos a crear x cuadrados desplazados 1 metros
for i in range(10):
    #print(i)

    #El ultimo node creado es
    last_node = ops.getNodeTags()[-1]
    #print(last_node)

    #Primero defino 4 nodos mas
    new_nodes = []
    for j, n in enumerate(nodes):
        ops.node(j+last_node+1, float(n[0]), float(n[1]), float(n[2]+1+i))
        new_nodes.append(j+last_node+1)

    #Ahora conecto los nodos nuevos entre si

    for i in new_nodes:

        j = i + 1
        
        if i == new_nodes[-1]:
            j = new_nodes[0]

            ops.element('Truss', element_id, i, j, A, 1, '-rho', gamma)
            element_id += 1
            elements.append((i,j))

        else:
            
            ops.element('Truss', element_id, i, j, A, 1, '-rho', gamma)
            element_id += 1
            elements.append((i,j))

    #Ahora conecto diagonalmente
    for i in old_nodes:

        if i == old_nodes[-1]:
            j = new_nodes[0]
            ops.element('Truss', element_id, i, j, A, 1, '-rho', gamma)
            element_id += 1
            elements.append((i,j))

        else:
            j = i + 5
            ops.element('Truss', element_id, i, j, A, 1, '-rho', gamma)
            element_id += 1
            elements.append((i,j))

        #print(i,j)

    #Ahora los conecto verticlamente
    for i in old_nodes:

        j = i + 4
        ops.element('Truss', element_id, i, j, A, 1)
        element_id += 1
        elements.append((i,j))
    
    #Ahora, los nodos antiguos los cambio antes de crear la nueva fila de nodos
    for i in range(len(new_nodes)):
        old_nodes[i] = new_nodes[i]



#Elements debe ser lista de listas y no tupla
elements_0_based = [(e[0], e[1]) for e in elements]
# Convertir elements_0_based a un array de NumPy
elements = np.array(elements_0_based)


#Genero el analisis
num_modes = 10
eigenvalues = ops.eigen(num_modes)

#Calculo las frecuencias naturales
eigenfrequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
print(f'Frecuencias naturales: {eigenfrequencies} Hz')

'''
#Comienzo con el grafico 3D
# Obtener las coordenadas de los nodos desde OpenSees
node_tags = ops.getNodeTags()
node_coords = [ops.nodeCoord(tag) for tag in node_tags]

# Crear el gráfico 3D con PyVista
plotter = pv.Plotter()
points = np.array(node_coords)

# Crear una nube de puntos para representar los nodos
point_cloud = pv.PolyData(points)
plotter.add_mesh(point_cloud, color='blue', point_size=10, render_points_as_spheres=True)

# Dibujar los elementos (barras) entre los nodos
for element in elements:
    start_node = np.array(ops.nodeCoord(element[0]))
    end_node = np.array(ops.nodeCoord(element[1]))
    line = pv.Line(start_node, end_node)
    plotter.add_mesh(line, color='red', line_width=3)

# Configurar y mostrar el gráfico
#plotter.show_grid()
#plotter.show()
'''


#Grafiquemos las frecuencias
def visualize_mode(mode_index, scale=0.2):

    global elements

    # Obtener todos los tags de nodos en el modelo
    node_tags = ops.getNodeTags()

    # Crear una lista para almacenar las coordenadas de cada nodo
    node_coords = []

    # Iterar sobre cada nodo y obtener sus coordenadas
    for tag in node_tags:
        coords = ops.nodeCoord(tag)  # Obtener las coordenadas (x, y, z) del nodo
        node_coords.append(coords)

    # Convertir la lista a un array de numpy
    nodes = np.array(node_coords)

    # Mostrar la variable 'nodes' creada
    plotter = pv.Plotter()

    # Get the mode shape for the specified mode
    mode_shape = np.array([ops.nodeEigenvector(n, mode_index) for n in ops.getNodeTags()])

    dl = np.sqrt(sum(((nodes[:, i].max() - nodes[:, i].min()) for i in range(3))))
    max_mode = abs(mode_shape[:, :]).max()

    scale_factor = scale * dl / max_mode

    # Apply the displacements (scaled) to the node coordinates
    displaced_nodes = nodes + scale_factor * mode_shape

    # Define connectivity (elements as 0-based indexing)

    elements = np.array(elements)  # Convierte la lista a un array de NumPy
    elements_0_based = elements - 1  # Restar 1 a todos los índices para que sean 0-based

    # Create the truss mesh (original configuration)
    truss = pv.PolyData(nodes)
    truss.lines = np.hstack([[2, e[0], e[1]] for e in elements_0_based])  # Set element connectivity
    truss["mode"] = "original"

    # Create the displaced mesh (deformed shape)
    truss_deformed = pv.PolyData(displaced_nodes)
    truss_deformed.lines = np.hstack([[2, e[0], e[1]] for e in elements_0_based])  # Set element connectivity
    truss_deformed["mode"] = f"Mode {mode_index}"

    # Plot original and deformed shapes
    plotter.add_mesh(truss, color='blue', label="Original")
    plotter.add_mesh(truss_deformed, color='red', label=f"Mode {mode_index}")

    # Add the frequency information to the plot
    frequency = eigenfrequencies[mode_index - 1]
    plotter.add_text(f"Mode {mode_index}: Frequency = {frequency:.4f} Hz", position='upper_left', font_size=10, color='black')

    plotter.add_legend()
    plotter.show()

# Visualize mode shapes
for i in range(1, num_modes + 1):
    visualize_mode(i, scale=0.1)  # Adjust scale to magnify displacements

#visualize_mode(1)

