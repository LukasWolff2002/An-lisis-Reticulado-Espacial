# -----------------------------------
# 1. Importación de Librerías
# -----------------------------------
import numpy as np
import pyvista as pv
import h5py
import sys
import os

# -----------------------------------
# 2. Definición de Funciones
# -----------------------------------

def cargar_modelo(file_path):
    """
    Carga el modelo de OpenSees desde un archivo HDF5.
    
    Parámetros:
        file_path (str): Ruta al archivo HDF5.
    
    Retorna:
        dict: Diccionario con los datos del modelo.
    """
    try:
        with h5py.File(file_path, 'r') as f:
            nodes_tags = f['nodes_tags'][()]
            nodes_xyz = f['nodes_xyz'][()]
            elements_tags = f['elements_tags'][()]
            elements_connectivities = f['elements_connectivities'][()]
        print(f"Modelo cargado exitosamente desde '{file_path}'.\n")
        return {
            'nodes_tags': nodes_tags.astype(int),
            'nodes_xyz': nodes_xyz.astype(float),
            'elements_tags': elements_tags.astype(int),
            'elements_connectivities': elements_connectivities.astype(int)
        }
    except Exception as e:
        print(f"Error al cargar el archivo HDF5: {e}")
        sys.exit(1)

def crear_malla_3d(modelo, escala=1.0, radio_cilindro=0.05, color_cilindro='blue', radio_nodo=0.1, color_nodo='red'):
    """
    Crea una malla 3D del modelo de OpenSees con un factor de escala.
    
    Parámetros:
        modelo (dict): Diccionario con datos del modelo.
        escala (float): Factor de escala para las coordenadas y radios.
        radio_cilindro (float): Radio de los cilindros que representan las barras (en metros).
        color_cilindro (str): Color de los cilindros.
        radio_nodo (float): Radio de las esferas que representan los nodos (en metros).
        color_nodo (str): Color de las esferas.
    
    Retorna:
        pyvista.PolyData: Malla 3D completa.
    """
    nodes = modelo['nodes_xyz'] * escala  # Aplicar escala a las coordenadas
    elements = modelo['elements_connectivities']
    
    # Crear una malla vacía
    malla = pv.PolyData()
    
    # Lista para almacenar todos los cilindros
    cilindros = []
    
    # Crear cilindros para cada barra
    for barra in elements:
        nodo_i, nodo_j = barra
        punto_i = nodes[nodo_i - 1]  # Asumiendo que los tags comienzan en 1
        punto_j = nodes[nodo_j - 1]
        
        # Crear una línea entre punto_i y punto_j
        linea = pv.Line(punto_i, punto_j)
        
        # Crear un cilindro (tubo) alrededor de la línea
        cilindro = linea.tube(radius=radio_cilindro * escala)  # Escalar el radio
        cilindro["color"] = color_cilindro
        cilindros.append(cilindro)
    
    # Combinar todos los cilindros en una sola malla
    if cilindros:
        malla = cilindros[0]
        for c in cilindros[1:]:
            malla = malla.merge(c)
    
    # Opcional: agregar esferas para representar nodos
    esferas = []
    for nodo, coord in zip(modelo['nodes_tags'], nodes):
        esfera = pv.Sphere(radius=radio_nodo * escala, center=coord)  # Escalar el radio
        esfera["color"] = color_nodo
        esferas.append(esfera)
    
    if esferas:
        malla = malla.merge(esferas[0])
        for e in esferas[1:]:
            malla = malla.merge(e)
    
    return malla

def visualizar_malla(malla):
    """
    Visualiza la malla 3D utilizando PyVista.
    
    Parámetros:
        malla (pyvista.PolyData): Malla 3D a visualizar.
    """
    if malla.n_points == 0:
        print("La malla está vacía. No hay nada para visualizar.")
        return
    
    plotter = pv.Plotter()
    plotter.add_mesh(malla, color='white', show_edges=False)
    plotter.add_axes()
    plotter.show()

def exportar_a_stl(malla, output_path):
    """
    Exporta la malla 3D a un archivo STL.
    
    Parámetros:
        malla (pyvista.PolyData): Malla 3D a exportar.
        output_path (str): Ruta donde se guardará el archivo STL.
    """
    try:
        malla.save(output_path)
        print(f"Malla exportada exitosamente a '{output_path}'.")
    except Exception as e:
        print(f"Error al exportar la malla a STL: {e}")

# -----------------------------------
# 3. Ejecución Principal
# -----------------------------------
def main():
    # Verificar si se proporcionó al menos un archivo como argumento
    if len(sys.argv) < 2:
        print("Uso: python modelo_3d.py <ruta_al_archivo_HDF5> [factor_de_escala]")
        sys.exit(1)
    
    file_path = sys.argv[1]
    
    # Verificar si se proporcionó un factor de escala
    if len(sys.argv) >= 3:
        try:
            escala = float(sys.argv[2])
            if escala <= 0:
                raise ValueError("El factor de escala debe ser un número positivo.")
        except ValueError as ve:
            print(f"Factor de escala inválido: {ve}")
            sys.exit(1)
    else:
        escala = 1.0  # Valor predeterminado
    
    print(f"Factor de escala: {escala}\n")
    
    # Cargar el modelo
    modelo = cargar_modelo(file_path)
    
    # Crear la malla 3D con el factor de escala
    malla = crear_malla_3d(modelo, escala=7, radio_cilindro=0.15, color_cilindro='blue', radio_nodo=0.1, color_nodo='red')
    
    # Visualizar la malla
    visualizar_malla(malla)
    
    # Definir la carpeta de salida y el nombre del archivo STL
    carpeta_destino = 'CODIGO/ENTREGA_3'
    nombre_stl = 'CajaSatelite.stl'
    ruta_stl = os.path.join(carpeta_destino, nombre_stl)
    
    # Crear la carpeta si no existe
    os.makedirs(carpeta_destino, exist_ok=True)
    
    # Exportar la malla a STL
    exportar_a_stl(malla, ruta_stl)

if __name__ == "__main__":
    main()
