import numpy as np
import pyvista as pv
import openseespy.opensees as ops
import h5py
#PARÁMETROS DE ENTRADA
limite1= 20 #Que tan largo desde la caja es la estructura
limite2= 36 #Que tan ancha es la estructura

gamma_fibra_carbono = 1.91 * 1000
E_fibra_carbono = 338e9
gamma_panel = 1.1
D1, D2 = 0.011 , 0.0001
A =np.pi*(D1**2-D2**2)/4

ops.wipe()  
ops.model('basic', '-ndm', 3, '-ndf', 3)

# Definición de materiales
ops.uniaxialMaterial('Elastic', 1, E_fibra_carbono)

def exportar_hf5():
    global nodos_totales, members, paneles_nodos
    barras=members
    # Necesito entregar un array con las coordenadas nodes_xyz
    # Otro array con los numeros de nodos nodes_tags
    # Otro array con los numeros de los elementos elements_tags
    # Otro array elements_section_info
    
    # Otro array elements_connectivities

    # Obtener los tags de los nodos y asegurarse de que sean enteros
    nodes = ops.getNodeTags()
    
    
    # Obtener las coordenadas de los nodos y asegurarse de que sean floats
    nodes_xyz = np.array([ops.nodeCoord(tag) for tag in nodes], dtype=float)
    nodes_tags = np.array(nodes, dtype=int)
    
    # Asegurarse de que nodos_fijos sea un array de enteros
    # Asumimos que nodos_fijos es una lista de listas: [node_tag, fix_x, fix_y, fix_z]
    nodos_fijos = []
    k=0
    for i, nodo in enumerate(nodes_tags):
        if i < 14:
            nodos_fijos.append([nodo, 1, 1, 1])  # Nodo fijo en x, y, z
        else:
            nodos_fijos.append([nodo, 0, 0, 0])  # Nodo libre en x, y, z
    
    # Convertir las fijaciones a un array de enteros
    nodes_fixities = np.array(nodos_fijos, dtype=int)
    
    # Crear un array con los tags de los elementos
    elements_tags = np.array([i + 1 for i in range(len(barras))], dtype=int)

    # Definir las conectividades de los elementos como enteros
    try:
        elements_connectivities = np.array([[int(barra[0]), int(barra[1])] for barra in barras], dtype=int)
    except ValueError as e:
        raise ValueError(f"Error al convertir conectividades de elementos: {e}")
    
    # Preparar la información de sección de los elementos
    elements_section_info = []
    for i in barras:
        elements_section_info.append([D2, (D1 - D2)/2])

    elements_section_info = np.array(elements_section_info, dtype=float)
    
    # Convertir panel_nodes a enteros
    panel_nodes = np.array(paneles_nodos , dtype=int)
    
    # Guardar en el archivo HDF5 con los tipos de datos correctos
    with h5py.File('CODIGO/ENTREGA_3/propuesta2.h5', 'w') as hf5:
        hf5.create_dataset('nodes_tags', data=nodes_tags)
        hf5.create_dataset('nodes_xyz', data=nodes_xyz)
        hf5.create_dataset('nodes_fixities', data=nodes_fixities)
        hf5.create_dataset('elements_tags', data=elements_tags)
        hf5.create_dataset('elements_connectivities', data=elements_connectivities)
        hf5.create_dataset('elements_section_info', data=elements_section_info)
        hf5.create_dataset('panel_nodes', data=panel_nodes)

    print("Datos exportados exitosamente a 'Caja_bernardo.h5'\n")

def plot_paneles_solares(plotter, nodos_paneles_solares1, nodos_paneles_solares2):

    # Procesar los paneles solares en ambas listas
    for i in range(2):  # Itera sobre las partes izquierda (0) y derecha (1) de cada lista
        lista1 = nodos_paneles_solares1[i]
        lista2 = nodos_paneles_solares2[i]

        for nodo in range(len(lista1) - 1):  # Itera hasta el penúltimo nodo para evitar errores de índice
            # Crear el primer triángulo intercalando los primeros nodos de lista1 y el primero de lista2
            if nodo < len(lista2):  # Asegura que lista2 tenga suficientes nodos
                nodos_panel1 = [lista1[nodo], lista1[nodo + 1], lista2[nodo]]
                
                surface1 = pv.PolyData(nodos_panel1)
                surface1.faces = [3, 0, 1, 2]  # Define la cara triangular
                plotter.add_mesh(surface1, color="gold", show_edges=True, opacity=0.7)

            # Crear el segundo triángulo intercalando los primeros nodos de lista2 y el siguiente de lista1
            if nodo + 1 < len(lista2):  # Asegura que no exceda el índice de lista2
                nodos_panel2 = [lista2[nodo], lista2[nodo + 1], lista1[nodo + 1]]
     
                surface2 = pv.PolyData(nodos_panel2)
                surface2.faces = [3, 0, 1, 2]  # Define la cara triangular
                plotter.add_mesh(surface2, color="gold", show_edges=True, opacity=0.7)
    
def masa_nodos(nodos_totales, members, gamma_fibra_carbono):
    masas_nodos = np.zeros(len(nodos_totales))  # Inicializamos una lista con ceros para almacenar la masa de cada nodo
    masa_estructura = 0  # Inicializamos la masa total de la estructura
    for member in members:
        nodo_i, nodo_j = member  # Extraer los nodos de cada miembro
        
        largo = np.linalg.norm(nodos_totales[nodo_i] - nodos_totales[nodo_j])
        
        masa_miembro = gamma_fibra_carbono * A * largo / 2
        
        masas_nodos[nodo_i] += masa_miembro
        masas_nodos[nodo_j] += masa_miembro
        
        masa_estructura += masa_miembro * 2

    print("Masa total de la estructura:", masa_estructura)
    return masas_nodos, masa_estructura

def calcular_masa_panel(paneles, nodos_totales):
    densidad_panel = 1.1  # Densidad del panel (kg/m^2)
    masas_nodos = np.zeros(len(nodos_totales))  # Inicializamos la masa de cada nodo en 0
    masa_total = 0  # Inicializamos la masa total de los paneles

    for panel in paneles:
        # Verifica que el panel tenga tres coordenadas
        if len(panel) == 3:
            # Obtiene las coordenadas de los tres puntos del panel
            punto1, punto2, punto3 = np.array(panel[0]), np.array(panel[1]), np.array(panel[2])

            # Calcula el área del triángulo en 3D
            vector1, vector2 = punto2 - punto1, punto3 - punto1
            area_panel = 0.5 * np.linalg.norm(np.cross(vector1, vector2))

            # Calcula la masa del panel
            masa_panel = densidad_panel * area_panel
            masa_total += masa_panel
            masa_nodo = masa_panel / 3  # Distribuye la masa equitativamente entre los tres nodos

            # Encuentra los índices de los nodos en nodos_totales y asigna la masa
            for punto in [punto1, punto2, punto3]:
                nodo_index = np.where((nodos_totales == punto).all(axis=1))[0]
                if nodo_index.size > 0:
                    masas_nodos[nodo_index[0]] += masa_nodo

        else:
            print(f"Advertencia: El panel {panel} no tiene exactamente tres puntos.")

    print("Masa total de los paneles solares:", masa_total)
    

    return masas_nodos, masa_total

def area_paneles_solares(nodos_paneles_solares1, nodos_paneles_solares2):
    paneles=[]
    # Procesar los paneles solares en ambas listas
    for i in range(2):  # Itera sobre las partes izquierda (0) y derecha (1) de cada lista
        lista1 = nodos_paneles_solares1[i]
        lista2 = nodos_paneles_solares2[i]

        for nodo in range(len(lista1) - 1):  # Itera hasta el penúltimo nodo para evitar errores de índice
            # Crear el primer triángulo intercalando los primeros nodos de lista1 y el primero de lista2
            if nodo < len(lista2):  # Asegura que lista2 tenga suficientes nodos
                nodos_panel1 = [lista1[nodo], lista1[nodo + 1], lista2[nodo]]
                paneles.append(nodos_panel1)

            # Crear el segundo triángulo intercalando los primeros nodos de lista2 y el siguiente de lista1
            if nodo + 1 < len(lista2):  # Asegura que no exceda el índice de lista2
                nodos_panel2 = [lista2[nodo], lista2[nodo + 1], lista1[nodo + 1]]
                paneles.append(nodos_panel2)

    area_total=0
    for panel in paneles:
        punto1 = np.array(panel[0])
        punto2 = np.array(panel[1])
        punto3 = np.array(panel[2])

        # Calcular el vector AB y AC
        vector_ab = punto2 - punto1
        vector_ac = punto3 - punto1

        # Calcular el producto cruzado de AB y AC
        producto_cruzado = np.cross(vector_ab, vector_ac)

        # Calcular el área como la mitad de la magnitud del producto cruzado
        area = 0.5 * np.linalg.norm(producto_cruzado)
        area_total+=area
    print("Área total de los paneles solares:", area_total)
    return paneles

#---------------------------------------------------------------------------------------------------------DEFINICIÓN DEL MODELO
# Definición inicial de nodos y apoyos
nodos_caja = np.array([
    [ 3.9, 3.3, 1.3],
    [-3.9, 3.3, 1.3],
    [-3.9, 3.3, -1.3],
    [ 3.9, 3.3,-1.3],
    [ 3.9,-3.3, 1.3],
    [-3.9,-3.3, 1.3],
    [-3.9,-3.3,-1.3],
    [ 3.9,-3.3,-1.3]
])
nodos_totales = nodos_caja
caras_caja = [
    [0, 3, 7, 4],   # Cara frontal
    [1, 2, 6, 5],   # Cara posterior
    [0, 1, 5, 4],   # Cara superior
    [2, 3, 7, 6],   # Cara inferior
    [0, 1, 2, 3],   # Cara izquierda
    [4, 5, 6, 7]    # Cara derecha
]

L0=1.3
L1=3.9
L2=3.3

apoyos_der = np.array([
    [ L1, L2,-L0],
    [-L1, L2,-L0],
    [   0,-L2,-L0],
    [ L1, L2, L0],
    [-L1, L2, L0],
    [   0,-L2, L0]
])

# Inicializar los nodos y miembros
nodos_totales = np.vstack((nodos_totales, apoyos_der))



# Definición de la función trusselator
def trusselatorder(limite1, limite2, nodos_totales, nodos_paneles_solares, L0, L1, L2):
    paso = 4
    m = len(nodos_totales)  # Contador de nodos, comienza donde terminan los nodos actuales
    for i in range(1, limite1 // paso+1):
        a = [
            [L1, L2, L0 + i * paso],
            [-L1, L2, L0 + i * paso],
            [0, -L2, L0 + i * paso] #añado el último
        ]
        if i==1:
            mbr = np.array([
                [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
                [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],  # Conecta con el triángulo anterior
                [m, m - 1], [m , m - 2], [m + 1, m - 1], [m+1,m-3], [m+2, m-3],[m+2,m-2]  # Conecta en diagonal
                 
            ])
            members =mbr
        else:
            mbr = np.array([
            [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
            [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],  # Conecta con el triángulo anterior
            [m, m - 1], [m , m - 2], [m + 1, m - 1], [m+1,m-3], [m+2, m-3],[m+2,m-2]  # Conecta en diagonal
        ])
        nodos_totales = np.vstack((nodos_totales, a))
        members = np.vstack((members, mbr))
        m += 3  # Incrementa el contador de nodos
    nodos_paneles_solares[0].append(a[2])
    nodos_paneles_solares[1].append(a[2])
    a=np.array([[0, 0, L0 + (i+1)*paso]])
    nodos_totales = np.vstack((nodos_totales, a))
    members = np.vstack((members, [[m, m-3], [m, m-2], [m, m-1]]))
    #----------------------------------------------------------------------------------------
    L0= L0 + (i)*paso
    m=m+1
    m_1=m-1
    for i in range(1, limite2//paso+1):
        a=[
                [L1+paso*i, L2, L0 ],
                [L1+paso*i, 0, L0+ paso],
                [L1 +paso*i, -L2, L0 ] #A este se le va a anclar los paneles solares
            ]
        nodos_paneles_solares[0].append(a[2])
        nodos_totales = np.vstack((nodos_totales, a))
        if i==1:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                            [m, m-1], [m, m-2], [m, m-3],
                            [m+1, m-1], [m+1, m-2],[m+1, m-4],
                            [m+2, m-1], [m+2, m-2], [m+2, m-4]])
                            #Se cierra el triángulo
        else:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                [m, m-3], [m+1, m-2], [m+2, m-1], #Se conecta con el triángulo anterior
                    [m, m-1], [m+1, m-3], [m+2, m-2]])
        m+=3
        members = np.vstack((members, mbr))
    
    for i in range(1, limite2//paso+1):
        a=[
                [-L1-paso*i, L2, L0 ],
                [-L1-paso*i, 0, L0+ paso],
                [-L1 -paso*i, -L2, L0 ] #A este se le va a anclar los paneles solares
            ]
        nodos_totales = np.vstack((nodos_totales, a))
        nodos_paneles_solares[1].append(a[2])
        if i==1:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                            [m, m_1], [m, m_1-1], [m, m_1-2],
                            [m+1, m_1], [m+1, m_1-1],[m+1, m_1-1],
                            [m+2, m_1], [m+2, m_1-1], [m+2, m_1-2]])
                            #Se cierra el triángulo
        else:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                [m, m-3], [m+1, m-2], [m+2, m-1], #Se conecta con el triángulo anterior
                    [m, m-1], [m+1, m-3], [m+2, m-2]])
        m+=3
        members = np.vstack((members, mbr))
   
    return nodos_totales, members, nodos_paneles_solares

def trusselatorizq(limite1, limite2, nodos_totales, members,nodos_paneles_solares, L0, L1, L2):
    paso = 4
    a1=8 # m-3
    b=9 # m-2
    c=10 # m-1
    m = len(nodos_totales)
    L0=-L0
    for i in range(1, limite1 // paso+1):
        a = [
            [L1, L2, L0 - i * paso],
            [-L1, L2, L0 - i * paso],
            [0, -L2, L0 - i * paso]
        ]
        if i==1:
            mbr = np.array([
                [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
                [m, a1], [m + 1, b], [m + 2, c],  # Conecta con el triángulo anterior
                [m, c], [m , b], [m + 1, c], [m+1,a1], [m+2, a1],[m+2,b]  # Conecta en diagonal
                 
            ])
        else:
            mbr = np.array([
            [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
            [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],  # Conecta con el triángulo anterior
            [m, m - 1], [m , m - 2], [m + 1, m - 1], [m+1,m-3], [m+2, m-3],[m+2,m-2]  # Conecta en diagonal
        ])
        nodos_totales = np.vstack((nodos_totales, a))
        members = np.vstack((members, mbr))
        m += 3  # Incrementa el contador de nodos
    nodos_paneles_solares[0].append(a[2])
    nodos_paneles_solares[1].append(a[2])
    a=np.array([[0, 0, L0 - (i+1)*paso]])
    nodos_totales = np.vstack((nodos_totales, a))
    members = np.vstack((members, [[m, m-3], [m, m-2], [m, m-1]]))
    #----------------------------------------------------------------------------------------
    L0= L0 - (i)*paso
    m=m+1
    m_1=m-1
    for i in range(1, limite2//paso+1):
        a=[
                [L1+paso*i, L2, L0 ],
                [L1+paso*i, 0, L0- paso],
                [L1 +paso*i, -L2, L0 ] #A este se le va a anclar los paneles solares
            ]
        nodos_paneles_solares[0].append(a[2])
        nodos_totales = np.vstack((nodos_totales, a))
        if i==1:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                            [m, m-1], [m, m-2], [m, m-3],
                            [m+1, m-1], [m+1, m-2],[m+1, m-4],
                            [m+2, m-1], [m+2, m-2], [m+2, m-4]])
                            #Se cierra el triángulo
        else:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                [m, m-3], [m+1, m-2], [m+2, m-1], #Se conecta con el triángulo anterior
                    [m, m-1], [m+1, m-3], [m+2, m-2]])
        m+=3
        members = np.vstack((members, mbr))
    
    for i in range(1, limite2//paso+1):
        a=[
                [-L1-paso*i, L2, L0 ],
                [-L1-paso*i, 0, L0- paso],
                [-L1 -paso*i, -L2, L0 ] #A este se le va a anclar los paneles solares
            ]
        nodos_paneles_solares[1].append(a[2])
        nodos_totales = np.vstack((nodos_totales, a))
        if i==1:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                            [m, m_1], [m, m_1-1], [m, m_1-2],
                            [m+1, m_1], [m+1, m_1-1],[m+1, m_1-1],
                            [m+2, m_1], [m+2, m_1-1], [m+2, m_1-2]])
                            #Se cierra el triángulo
        else:
            mbr= np.array([[m, m+1], [m+1, m+2], [m+2, m], #Se cierra el triángulo
                [m, m-3], [m+1, m-2], [m+2, m-1], #Se conecta con el triángulo anterior
                    [m, m-1], [m+1, m-3], [m+2, m-2]])
        m+=3
        members = np.vstack((members, mbr))
   
    return nodos_totales, members, nodos_paneles_solares

# Llamar a la función trusselator
nodos_paneles_solares=[[],[]]
nodos_totales, members, nodos_paneles_solares1 = trusselatorder(limite1, limite2, nodos_totales, nodos_paneles_solares, L0, L1, L2)


nodos_paneles_solares=[[],[]]
nodos_totales, members, nodos_paneles_solares2 = trusselatorizq(limite1, limite2, nodos_totales, members, nodos_paneles_solares, L0, L1, L2)


nodos_paneles=[]
for i in range(len(nodos_paneles_solares1[0])):
    nodos_paneles.append(nodos_paneles_solares1[0][i])

#---------------------------------------------------------------------------------------------------------DEFINICIÓN DEL MODELO
#---------------------------------------------------------------------------------------------------------EJECUCIÓN DEL MODELO
for i, nodo in enumerate(nodos_totales):
    ops.node(i + 1, *nodo)

for i, member in enumerate(members):
    nodo_i = int(member[0] + 1)
    nodo_j = int(member[1] + 1)
    ops.element('Truss', i + 1, nodo_i, nodo_j, A, 1, '-rho', gamma_fibra_carbono)



# Fijar los nodos de los apoyos en todos los grados de libertad (x, y, z)
for idx in range(14):
    ops.fix(idx + 1, 1, 1, 1)  # Fijar todos los grados de libertad en el apoyo derecho


#---------------------------------------------------------------------------------------------------------Masa de los nodos
# Calcular la masa de los nodos de la estructura
masas_nodos, masa_estructura = masa_nodos(nodos_totales, members, gamma_fibra_carbono)

# Calcular la masa de los paneles solares
paneles = area_paneles_solares(nodos_paneles_solares1, nodos_paneles_solares2)
masas_nodo_panel, masa_total_paneles = calcular_masa_panel(paneles, nodos_totales)

# Sumar la masa de los paneles a la masa estructural de cada nodo

# Asignar las masas combinadas en OpenSees
for i, masa in enumerate(masas_nodos, start=1):
    ops.mass(i, masa, masa, masa)


RME= masa_estructura/masa_total_paneles

print("RME de la estructura:", RME)

num_modos = 3  # Número de modos a calcular
eigenvalues = ops.eigen(num_modos)  # Realiza el análisis modal
frecuencias = [np.sqrt(eigenval) / (2 * np.pi) for eigenval in eigenvalues]
for i, freq in enumerate(frecuencias, start=1):
    print(f"Frecuencia del modo {i}: {freq:.8f} Hz")

paneles_nodos = []

# Obtener todos los nodos del modelo
nodes = ops.getNodeTags()
nodes_coords = {node: np.array(ops.nodeCoord(node)) for node in nodes}  # Diccionario de nodos y sus coordenadas

# Recorrer los paneles y encontrar los nodos correspondientes
tolerance = 1e-6  # Tolerancia para comparación
for panel in paneles:
    panel_nodos = []
    for coord in panel:  # Coordenadas de los puntos del panel
        for node, node_coords in nodes_coords.items():  # Comparar con las coordenadas de los nodos
            if np.linalg.norm(node_coords - np.array(coord)) < tolerance:
                panel_nodos.append(node)  # Añadir el nodo correspondiente
                break
    paneles_nodos.append(panel_nodos)

#---------------------------------------------------------------------------------------------------------EJECUCIÓN DEL MODELO
plotter = pv.Plotter()

# Agregar nodos como puntos en el gráfico
plotter.add_points(nodos_totales, color="black", point_size=10, render_points_as_spheres=True)

# Agregar las caras de la caja como superficies
for cara in caras_caja:
    puntos_panel = nodos_caja[cara]
    surface = pv.PolyData(puntos_panel)
    surface.faces = [4, 0, 1, 2, 3]
    plotter.add_mesh(surface, color='lightblue', show_edges=True, opacity=0.5)

# Agregar los elementos como líneas
for member in members:
    start_node = nodos_totales[member[0]]
    end_node = nodos_totales[member[1]]
    line = pv.Line(start_node, end_node)
    plotter.add_mesh(line, color="red", line_width=2)

plot_paneles_solares(plotter, nodos_paneles_solares1, nodos_paneles_solares2)

# Agregar etiquetas a todos los nodos
for idx, nodo in enumerate(nodos_totales):
    plotter.add_point_labels([nodo], [f"{idx+1}"], font_size=10, text_color="blue")



# Visualización final
plotter.add_axes()
plotter.show()





exportar_hf5()
