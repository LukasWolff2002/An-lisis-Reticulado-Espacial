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
    nodos_fijos_int = []
    for fij in nodos_totales:
        # Convertir cada elemento a entero
        try:
            fij_int = [int(f) for f in fij]
            nodos_fijos_int.append(fij_int)
        except ValueError as e:
            raise ValueError(f"Error al convertir restricciones de nodo {fij}: {e}")
    
    nodes_fixities = np.array(nodos_fijos_int, dtype=int)
    
    # Definir los tags de los elementos como enteros
    elements_tags = np.array([i for i in range(1, len(barras)+1)], dtype=int)
    
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
    with h5py.File('Caja_bernardo.h5', 'w') as hf5:
        hf5.create_dataset('nodes_tags', data=nodes_tags)
        hf5.create_dataset('nodes_xyz', data=nodes_xyz)
        hf5.create_dataset('nodes_fixities', data=nodes_fixities)
        hf5.create_dataset('elements_tags', data=elements_tags)
        hf5.create_dataset('elements_connectivities', data=elements_connectivities)
        hf5.create_dataset('elements_section_info', data=elements_section_info)
        hf5.create_dataset('panel_nodes', data=panel_nodes)

    print("Datos exportados exitosamente a 'Caja_bernardo.h5'\n")