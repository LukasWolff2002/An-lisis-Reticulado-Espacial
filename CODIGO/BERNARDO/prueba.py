import numpy as np
import pyvista as pv

# Definición inicial de nodos y apoyos
nodos_caja = np.array([
    [3.3, 3.9, 1.3],
    [3.3, -3.9, 1.3],
    [3.3, -3.9, -1.3],
    [3.3, 3.9, -1.3],
    [-3.3, 3.9, 1.3],
    [-3.3, -3.9, 1.3],
    [-3.3, -3.9, -1.3],
    [3.3, 3.9, -1.3]
])

apoyos_der = np.array([
    [3.2, 3.8, 1.3],
    [-3.2, 3.8, 1.3],
    [0, -3.8, 1.3]
])

# Inicializar los nodos y miembros
nodos_totales = np.vstack((nodos_caja, apoyos_der))
members = np.array([[9, 10], [10, 11], [11, 9]])

# Definición de la función trusselatorder (lado derecho)
def trusselatorder(limite1, limite2, nodos_totales, members):
    paso = 4
    m = len(nodos_totales)  # Contador de nodos
    for i in range(1, limite1 // paso + 1):
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
    a = np.array([[0, 0, L0 + (i + 1) * paso]])
    nodos_totales = np.vstack((nodos_totales, a))
    members = np.vstack((members, [[m, m - 3], [m, m - 2], [m, m - 1]]))
    return nodos_totales, members

# Definición de la función trusselatorizq (lado izquierdo)
def trusselatorizq(limite1, limite2, nodos_totales, members):
    paso = 4
    m = len(nodos_totales)  # Contador de nodos
    for i in range(1, limite1 // paso + 1):
        L0 = 1.3
        L1 = 3.2
        L2 = 3.8
        a = [
            [L1, L2, -L0 - i * paso],
            [-L1, L2, -L0 - i * paso],
            [0, -L2, -L0 - i * paso]
        ]
        mbr = np.array([
            [m, m + 1], [m + 1, m + 2], [m + 2, m],  # Cierra el triángulo
            [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],  # Conecta con el triángulo anterior
            [m, m - 1], [m + 1, m - 3], [m + 2, m - 2]  # Conecta en diagonal
        ])
        nodos_totales = np.vstack((nodos_totales, a))
        members = np.vstack((members, mbr))
        m += 3  # Incrementa el contador de nodos
    a = np.array([[0, 0, -L0 - (i + 1) * paso]])
    nodos_totales = np.vstack((nodos_totales, a))
    members = np.vstack((members, [[m, m - 3], [m, m - 2], [m, m - 1]]))
    return nodos_totales, members

# Generar nodos y miembros para ambos lados
nodos_totales, members = trusselatorder(40, 40, nodos_totales, members)
nodos_totales, members = trusselatorizq(40, 40, nodos_totales, members)

# Ajuste de índices: restar 1 para compatibilidad con NumPy
members -= 1

# Visualización con PyVista
plotter = pv.Plotter()

# Agregar nodos como puntos en el gráfico
plotter.add_points(nodos_totales, color="black", point_size=10, render_points_as_spheres=True)

# Agregar los elementos como líneas
for member in members:
    start_node = nodos_totales[member[0]]
    end_node = nodos_totales[member[1]]
    line = pv.Line(start_node, end_node)
    plotter.add_mesh(line, color="red", line_width=2)

# Visualización final
plotter.add_axes()
plotter.show()
