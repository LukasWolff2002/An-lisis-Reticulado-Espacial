import matplotlib.pyplot as plt

# Datos de las tablas
tables = {
    "grafico_frecuencias_2vanos": {
        "Frecuencia (Hz)": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "A_10%": [0.245, 0.283, 0.314, 0.666, 0.673, 0.744, 1.014, 1.046, 1.206, 1.454],
        "A_20%": [0.349, 0.400, 0.445, 0.947, 0.958, 1.054, 1.442, 1.488, 1.717, 2.070],
        "A_30%": [0.425, 0.484, 0.540, 1.154, 1.166, 1.280, 1.757, 1.812, 2.092, 2.522]
    },
    "grafico_frecuencias_4vanos": {
        "Frecuencia (Hz)": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "A_10%": [0.103, 0.141, 0.162, 0.342, 0.393, 0.449, 0.597, 0.624, 0.688, 0.735],
        "A_20%": [0.146, 0.195, 0.228, 0.486, 0.553, 0.633, 0.849, 0.887, 0.974, 1.039],
        "A_30%": [0.175, 0.228, 0.272, 0.582, 0.656, 0.754, 1.015, 1.060, 1.162, 1.236]
    },
    "grafico_frecuencias_8vanos": {
        "Frecuencia (Hz)": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "A_10%": [0.035, 0.064, 0.081, 0.147, 0.195, 0.234, 0.307, 0.336, 0.388, 0.444],
        "A_20%": [0.049, 0.083, 0.115, 0.208, 0.266, 0.330, 0.435, 0.474, 0.548, 0.629],
        "A_30%": [0.059, 0.095, 0.139, 0.253, 0.314, 0.399, 0.528, 0.574, 0.663, 0.760]
    },
    "grafico_frecuencias_16vanos": {
        "Frecuencia (Hz)": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "A_10%": [0.010, 0.025, 0.040, 0.052, 0.086, 0.113, 0.130, 0.161, 0.192, 0.212],
        "A_20%": [0.014, 0.029, 0.057, 0.074, 0.112, 0.160, 0.184, 0.222, 0.270, 0.299],
        "A_30%": [0.017, 0.031, 0.069, 0.090, 0.130, 0.194, 0.223, 0.264, 0.326, 0.363]
    }
}

# Generar y guardar gráficos
for filename, data in tables.items():
    plt.figure()
    plt.plot(data["Frecuencia (Hz)"], data["A_10%"], label="$A_{10\\%}$", marker='o')
    plt.plot(data["Frecuencia (Hz)"], data["A_20%"], label="$A_{20\\%}$", marker='o')
    plt.plot(data["Frecuencia (Hz)"], data["A_30%"], label="$A_{30\\%}$", marker='o')
    
    # Etiquetas y título
    plt.xlabel("Modo")
    plt.ylabel("Frecuencia (Hz)")
    plt.title(filename.replace("_", " ").capitalize())
    plt.legend()
    
    # Guardar el gráfico
    plt.savefig(f"C:/Users/pedro/OneDrive/Escritorio/Semestre VIII/MCOC/p2/p2e1/Analisis-Reticulado-Espacial/{filename}.png")
    plt.close()
