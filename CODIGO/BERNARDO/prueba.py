import numpy as np
import pyvista as pv
import openseespy.opensees as ops

# Parámetros de entrada
limite1 = 20  # Largo desde la caja de la estructura
limite2 = 36  # Ancho de la estructura

gamma_fibra_carbono = 1.91 * 1000  # Densidad de la fibra de carbono (kg/m^3)
E_fibra_carbono = 338e9            # Módulo de elasticidad de la fibra de carbono (Pa)
gamma_panel = 1.1                  # Densidad de los paneles solares (kg/m^2)
D1, D2 = 0.10, 0.006               # Diámetros exterior e interior (m)
A = np.pi * (D1**2 - D2**2) / 4    # Área de la sección transversal

# Inicialización del modelo en OpenSeesPy
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# Definición de materiales
ops.uniaxialMaterial('Elastic', 1, E_fibra_carbono)

def plot_paneles_solares(plotter, nodos_paneles_solares1, nodos_paneles_solares2):
    # Procesar los paneles solares en ambas listas
    for i in range(2):  # Itera sobre las partes izquierda (0) y derecha (1)
        lista1 = nodos_paneles_solares1[i]
        lista2 = nodos_paneles_solares2[i]

        for nodo in range(len(lista1) - 1):
            # Crear el primer triángulo
            if nodo < len(lista2):
                nodos_panel1 = [lista1[nodo], lista1[nodo + 1], lista2[nodo]]
                surface1 = pv.PolyData(nodos_panel1)
                surface1.faces = [3, 0, 1, 2]
                plotter.add_mesh(surface1, color="gold", show_edges=True, opacity=0.7)

            # Crear el segundo triángulo
            if nodo + 1 < len(lista2):
                nodos_panel2 = [lista2[nodo], lista2[nodo + 1], lista1[nodo + 1]]
                surface2 = pv.PolyData(nodos_panel2)
                surface2.faces = [3, 0, 1, 2]
                plotter.add_mesh(surface2, color="gold", show_edges=True, opacity=0.7)

def masa_nodos(nodos_totales, members, gamma_fibra_carbono):
    masas_nodos = np.zeros(len(nodos_totales))
    A = np.pi * (D1**2 - D2**2) / 4
    masa_estructura = 0
    for member in members:
        nodo_i, nodo_j = member
        largo = np.linalg.norm(nodos_totales[nodo_i] - nodos_totales[nodo_j])
        masa_miembro = gamma_fibra_carbono * A * largo / 2
        masas_nodos[nodo_i] += masa_miembro
        masas_nodos[nodo_j] += masa_miembro
        masa_estructura += masa_miembro * 2

    print("Masa total de la estructura:", masa_estructura)
    return masas_nodos, masa_estructura

def calcular_masa_panel(paneles, nodos_totales):
    densidad_panel = gamma_panel  # Densidad del panel (kg/m^2)
    masas_nodos = np.zeros(len(nodos_totales))
    masa_total = 0

    for panel in paneles:
        if len(panel) == 3:
            punto1, punto2, punto3 = np.array(panel[0]), np.array(panel[1]), np.array(panel[2])
            vector1, vector2 = punto2 - punto1, punto3 - punto1
            area_panel = 0.5 * np.linalg.norm(np.cross(vector1, vector2))
            masa_panel = densidad_panel * area_panel
            masa_total += masa_panel
            masa_nodo = masa_panel / 3

            # Encuentra los índices de los nodos
            for punto in [punto1, punto2, punto3]:
                nodo_index = np.where((nodos_totales == punto).all(axis=1))[0]
                if nodo_index.size > 0:
                    masas_nodos[nodo_index[0]] += masa_nodo
        else:
            print(f"Advertencia: El panel {panel} no tiene exactamente tres puntos.")

    print("Masa total de los paneles solares:", masa_total)
    return masas_nodos, masa_total

def area_paneles_solares(nodos_paneles_solares1, nodos_paneles_solares2):
    paneles = []
    # Procesar los paneles solares en ambas listas
    for i in range(2):  # Partes izquierda y derecha
        lista1 = nodos_paneles_solares1[i]
        lista2 = nodos_paneles_solares2[i]

        for nodo in range(len(lista1) - 1):
            if nodo < len(lista2):
                nodos_panel1 = [lista1[nodo], lista1[nodo + 1], lista2[nodo]]
                paneles.append(nodos_panel1)
            if nodo + 1 < len(lista2):
                nodos_panel2 = [lista2[nodo], lista2[nodo + 1], lista1[nodo + 1]]
                paneles.append(nodos_panel2)

    area_total = 0
    for panel in paneles:
        punto1, punto2, punto3 = np.array(panel[0]), np.array(panel[1]), np.array(panel[2])
        vector_ab = punto2 - punto1
        vector_ac = punto3 - punto1
        producto_cruzado = np.cross(vector_ab, vector_ac)
        area = 0.5 * np.linalg.norm(producto_cruzado)
        area_total += area

    print("Área total de los paneles solares:", area_total)
    return paneles

# Definición de los nodos de la caja central
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

# Nodos de los apoyos derecho e izquierdo
apoyos_der = np.array([
    [3.2, 3.8, 1.3],
    [-3.2, 3.8, 1.3],
    [0, -3.8, 1.3]
])
apoyos_izq = np.array([
    [3.2, 3.8, -1.3],
    [-3.2, 3.8, -1.3],
    [0, -3.8, -1.3]
])

# Inicializar los nodos y miembros
nodos_totales = np.vstack((nodos_caja, apoyos_der))
members = np.array([[8, 9], [9, 10], [10, 8]])  # Ajustado para índices correctos

# Función para generar la estructura de celosía derecha
def trusselator(limite1, limite2, nodos_totales, members, nodos_paneles_solares, direccion='derecha'):
    paso = 4
    L0 = 1.3 if direccion == 'derecha' else -1.3
    signo = 1 if direccion == 'derecha' else -1
    L1 = 3.2
    L2 = 3.8
    m = len(nodos_totales)

    # Generación en la dirección longitudinal
    for i in range(1, limite1 // paso + 1):
        a = [
            [L1, L2, L0 + signo * i * paso],
            [-L1, L2, L0 + signo * i * paso],
            [0, -L2, L0 + signo * i * paso]
        ]
        mbr = np.array([
            [m, m + 1], [m + 1, m + 2], [m + 2, m],
            [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],
            [m, m - 1], [m + 1, m - 3], [m + 2, m - 2]
        ])
        nodos_totales = np.vstack((nodos_totales, a))
        members = np.vstack((members, mbr))
        m += 3

    # Nodo central
    a = np.array([[0, 0, L0 + signo * (i + 1) * paso]])
    nodos_totales = np.vstack((nodos_totales, a))
    members = np.vstack((members, [[m, m - 3], [m, m - 2], [m, m - 1]]))
    m += 1

    # Generación en la dirección transversal
    for i in range(1, limite2 // paso + 1):
        a = [
            [L1 + paso * i, L2, L0 + signo * i * paso],
            [L1 + paso * i, 0, L0 + signo * (i + 1) * paso],
            [L1 + paso * i, -L2, L0 + signo * i * paso]
        ]
        nodos_paneles_solares[0].append(a[2])
        nodos_totales = np.vstack((nodos_totales, a))
        if i == 1:
            mbr = np.array([
                [m, m + 1], [m + 1, m + 2], [m + 2, m],
                [m, m - 1], [m, m - 2], [m, m - 3],
                [m + 1, m - 1], [m + 1, m - 2], [m + 1, m - 4],
                [m + 2, m - 1], [m + 2, m - 2], [m + 2, m - 4]
            ])
        else:
            mbr = np.array([
                [m, m + 1], [m + 1, m + 2], [m + 2, m],
                [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],
                [m, m - 1], [m + 1, m - 3], [m + 2, m - 2]
            ])
        members = np.vstack((members, mbr))
        m += 3

    for i in range(1, limite2 // paso + 1):
        a = [
            [-L1 - paso * i, L2, L0 + signo * i * paso],
            [-L1 - paso * i, 0, L0 + signo * (i + 1) * paso],
            [-L1 - paso * i, -L2, L0 + signo * i * paso]
        ]
        nodos_paneles_solares[1].append(a[2])
        nodos_totales = np.vstack((nodos_totales, a))
        if i == 1:
            mbr = np.array([
                [m, m + 1], [m + 1, m + 2], [m + 2, m],
                [m, m - 1], [m, m - 2], [m, m - 3],
                [m + 1, m - 1], [m + 1, m - 2], [m + 1, m - 4],
                [m + 2, m - 1], [m + 2, m - 2], [m + 2, m - 4]
            ])
        else:
            mbr = np.array([
                [m, m + 1], [m + 1, m + 2], [m + 2, m],
                [m, m - 3], [m + 1, m - 2], [m + 2, m - 1],
                [m, m - 1], [m + 1, m - 3], [m + 2, m - 2]
            ])
        members = np.vstack((members, mbr))
        m += 3

    return nodos_totales, members, nodos_paneles_solares

# Generación de la celosía derecha
nodos_paneles_solares1 = [[], []]
nodos_totales, members, nodos_paneles_solares1 = trusselator(
    limite1, limite2, nodos_totales, members, nodos_paneles_solares1, direccion='derecha'
)

# Añadir los apoyos izquierdos
nodos_totales = np.vstack((nodos_totales, apoyos_izq))
indices_apoyos_izq = [len(nodos_totales) - 3, len(nodos_totales) - 2, len(nodos_totales) - 1]
members = np.vstack((members, [[len(nodos_totales) - 3, len(nodos_totales) - 2],
                               [len(nodos_totales) - 2, len(nodos_totales) - 1],
                               [len(nodos_totales) - 1, len(nodos_totales) - 3]]))

# Generación de la celosía izquierda
nodos_paneles_solares2 = [[], []]
nodos_totales, members, nodos_paneles_solares2 = trusselator(
    limite1, limite2, nodos_totales, members, nodos_paneles_solares2, direccion='izquierda'
)

# Definir todos los nodos en OpenSeesPy
for i, nodo in enumerate(nodos_totales):
    ops.node(i + 1, *nodo)

# Definir los elementos truss
for i, member in enumerate(members):
    nodo_i = int(member[0] + 1)
    nodo_j = int(member[1] + 1)
    ops.element('Truss', i + 1, nodo_i, nodo_j, A, 1, '-rho', gamma_fibra_carbono)

# Índices de los nodos de apoyo derecho
indices_apoyos_der = [8, 9, 10]

# Fijar los nodos de los apoyos
for idx in indices_apoyos_der:
    ops.fix(idx + 1, 1, 1, 1)

for idx in indices_apoyos_izq:
    ops.fix(idx + 1, 1, 1, 1)

# Cálculo de masas nodales debido a los elementos estructurales
masas_nodos_estructura, masa_estructura = masa_nodos(nodos_totales, members, gamma_fibra_carbono)

# Asignación de masas de la estructura a los nodos
for i, masa in enumerate(masas_nodos_estructura, start=1):
    ops.mass(i, masa, masa, masa)

# Cálculo de áreas y masas de los paneles solares
paneles = area_paneles_solares(nodos_paneles_solares1, nodos_paneles_solares2)
masas_nodos_paneles, masa_total_paneles = calcular_masa_panel(paneles, nodos_totales)

# Sumar las masas de los paneles a las masas existentes en los nodos
for i, masa in enumerate(masas_nodos_paneles, start=1):
    masa_existente = ops.nodeMass(i)
    masa_total_nodo = masa_existente[0] + masa
    ops.mass(i, masa_total_nodo, masa_total_nodo, masa_total_nodo)

# Cálculo de la relación masa estructura / masa paneles
RME = masa_estructura / masa_total_paneles
print("RME de la estructura:", RME)

# Análisis modal
eigenvalues = ops.eigen(3)

# Conversión de valores propios a frecuencias naturales
frecuencias = [np.sqrt(eigenval) / (2 * np.pi) for eigenval in eigenvalues]

# Mostrar las frecuencias de los primeros tres modos de vibración
for i, freq in enumerate(frecuencias, start=1):
    print(f"Frecuencia del modo {i}: {freq:.4f} Hz")

# Visualización con PyVista
plotter = pv.Plotter()

# Añadir nodos al gráfico
plotter.add_points(nodos_totales, color="black", point_size=10, render_points_as_spheres=True)

# Dibujar las caras de la caja central
for cara in caras_caja:
    puntos_panel = nodos_caja[cara]
    surface = pv.PolyData(puntos_panel)
    surface.faces = [4, 0, 1, 2, 3]
    plotter.add_mesh(surface, color='lightblue', show_edges=True, opacity=0.5)

# Dibujar los elementos estructurales
for member in members:
    start_node = nodos_totales[member[0]]
    end_node = nodos_totales[member[1]]
    line = pv.Line(start_node, end_node)
    plotter.add_mesh(line, color="red", line_width=2)

# Dibujar los paneles solares
plot_paneles_solares(plotter, nodos_paneles_solares1, nodos_paneles_solares2)

# Etiquetar los nodos de los apoyos
for idx in indices_apoyos_der:
    plotter.add_point_labels([nodos_totales[idx]], [f"Apoyo Der {idx+1}"], font_size=15, text_color="green")

for idx in indices_apoyos_izq:
    plotter.add_point_labels([nodos_totales[idx]], [f"Apoyo Izq {idx+1}"], font_size=15, text_color="green")

# Mostrar el gráfico
plotter.add_axes()
plotter.show()
