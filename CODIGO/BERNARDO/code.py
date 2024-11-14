import numpy as np
import pyvista as pv
import openseespy.opensees as ops
#PARÁMETROS DE ENTRADA
limite1= 20 #Que tan largo desde la caja es la estructura
limite2= 48 #Que tan ancha es la estructura

gamma_fibra_carbono = 1.91 * 1000
E_fibra_carbono = 338e9
gamma_panel = 1.1
D1, D2 = 0.010, 0.006
A =np.pi*(D1**2-D2**2)/4

ops.wipe()  
ops.model('basic', '-ndm', 3, '-ndf', 3)

# Definición de materiales
ops.uniaxialMaterial('Elastic', 1, E_fibra_carbono)

#---------------------------------------------------------------------------------------------------------DEFINICIÓN DEL MODELO
# Definición inicial de nodos y apoyos
nodos_caja = np.array([
    [3.3, 3.9, 1.3],
    [3.3, -3.9, 1.3],
    [3.3, -3.9, -1.3],
    [3.3, 3.9, -1.3],
    [-3.3, 3.9, 1.3],
    [-3.3, -3.9, 1.3],
    [-3.3, -3.9, -1.3],
    [-3.3, 3.9, -1.3]
])
caras_caja = [
    [0, 3, 7, 4],   # Cara frontal
    [1, 2, 6, 5],   # Cara posterior
    [0, 1, 5, 4],   # Cara superior
    [2, 3, 7, 6],   # Cara inferior
    [0, 1, 2, 3],   # Cara izquierda
    [4, 5, 6, 7]    # Cara derecha
]

apoyos_der = np.array([
    [3.2, 3.8, 1.3],
    [-3.2, 3.8, 1.3],
    [0, -3.8, 1.3]
])

# Inicializar los nodos y miembros
nodos_totales = np.vstack((nodos_caja, apoyos_der))
members = np.array([[9, 10], [10, 11], [11, 9]])

# Definición de la función trusselator
def trusselatorder(limite1, limite2, nodos_totales, members):
    paso = 4
    m = len(nodos_totales)  # Contador de nodos, comienza donde terminan los nodos actuales
    for i in range(1, limite1 // paso+1):
        L0 = 1.3
        L1 = 3.2
        L2 = 3.8
        a = [
            [L1, L2, L0 + i * paso],
            [-L1, L2, L0 + i * paso],
            [0, -L2, L0 + i * paso]
        ]
        mbr = np.array([
            [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
            [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],  # Conecta con el triángulo anterior
            [m, m - 1], [m + 1, m - 3], [m + 2, m - 2]  # Conecta en diagonal
        ])
        nodos_totales = np.vstack((nodos_totales, a))
        members = np.vstack((members, mbr))
        m += 3  # Incrementa el contador de nodos
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
                [L1 +paso*i, -L2, L0 ]
            ]
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
                [-L1 -paso*i, -L2, L0 ]
            ]
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
   
    return nodos_totales, members

def trusselatorizq(limite1, limite2, nodos_totales, members):
    paso = 4
    m = len(nodos_totales) 
    for i in range(1, limite1 // paso+1):
        L0 = -1.3
        L1 = 3.2
        L2 = 3.8
        a = [
            [L1, L2, L0 - i * paso],
            [-L1, L2, L0 - i * paso],
            [0, -L2, L0 - i * paso]
        ]
        mbr = np.array([
            [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
            [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],  # Conecta con el triángulo anterior
            [m, m - 1], [m + 1, m - 3], [m + 2, m - 2]  # Conecta en diagonal
        ])
        nodos_totales = np.vstack((nodos_totales, a))
        members = np.vstack((members, mbr))
        m += 3  # Incrementa el contador de nodos
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
                [L1 +paso*i, -L2, L0 ]
            ]
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
                [-L1 -paso*i, -L2, L0 ]
            ]
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
   
    return nodos_totales, members

# Llamar a la función trusselator

nodos_totales, members = trusselatorder(limite1, limite2, nodos_totales, members)

m= len(nodos_totales)+1

apoyos_der = np.array([
    [3.2, 3.8, -1.3],
    [-3.2, 3.8, -1.3],
    [0, -3.8, -1.3]
])

mbr= [[m, m+1], [m+1, m+2], [m+2, m]] 

nodos_totales=np.vstack((nodos_totales, apoyos_der))
members=np.vstack((members, mbr))

nodos_totales, members = trusselatorizq(limite1, limite2, nodos_totales, members)
#---------------------------------------------------------------------------------------------------------DEFINICIÓN DEL MODELO
#---------------------------------------------------------------------------------------------------------EJECUCIÓN DEL MODELO
for i, nodo in enumerate(nodos_totales):
    ops.node(i + 1, *nodo)

for i, member in enumerate(members):
    nodo_i = int(member[0] + 1)
    nodo_j = int(member[1] + 1)
    ops.element('Truss', i + 1, nodo_i, nodo_j, A, 1, '-rho', gamma_fibra_carbono)
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

# Visualización final
plotter.add_axes()
plotter.show()

