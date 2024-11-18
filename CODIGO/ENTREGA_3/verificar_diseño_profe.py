#!/bin//bash
import sys
sys.path.append("/home/jaabell/Repositories/OpenSees-ASDPlastic/build")
import openseespy.opensees as ops  #es distinto para ustedesimport numpy as np
import numpy as np
import h5py

# Condiciones de diseño
# ===================================================================================

g = 9.81                    # m/s²

# Potencia 1MW (considerar 300 W/m² de eficiencia de panel solar)
eficiencia_panel = 300      # W /m²
potencia_requerida = 1e6    # W
area_min = potencia_requerida / eficiencia_panel

# Peso panel solar: 1.1 kg / m²
ρ_panel = 1.1           #kg / m²

# Fracción de masa < 20%
RME_max = 0.2

# Frecuencia mínima del primer modo > 0.1Hz
f_min_requerida = 0.1

# Aceleración máxima 0.1g (en cualquier dirección)
a_max = 0.1 * g 

# Factor de seguridad = 2.0
FS = 2.0

# Desangulación < 2° para un cambio de temperatura de 150°C en la zona del panel solar.
θ_max = 2
ΔT = 150.               #°C 

# Material compuesto en base a fibra de carbono de alto módulo M55J (especificaciones Download especificaciones)
alpha = 1e-6             # CTE M55J - 1/°C
E_compuesto = 338e9      # Pa
σ_max_compuesto = 2020e6 # Pa
ρ_compuesto = 1200       # kg /m³

# Tamaño satélite: caja de 6.6 m x 2.6 m x 7.8 m
# El reticulado se puede apoyar de cualquier forma en esta "caja".
Lx = 7.8
Ly = 6.6
Lz = 2.6

if len(sys.argv) == 1:
    fname = "CODIGO/ENTREGA_3/propuesta1.h5"
else:
    fname = sys.argv[1]

print(f"Revisando el diseño en {fname}")
print("")

#Leer modelo
f = h5py.File(fname)

nodes_tags = f["nodes_tags"][()]
nodes_xyz = f["nodes_xyz"][()]
nodes_fixities = f["nodes_fixities"][()]
elements_tags = f["elements_tags"][()]
elements_connectivities = f["elements_connectivities"][()]
elements_section_info = f["elements_section_info"][()]
panel_nodes = f["panel_nodes"][()]

f.close()


#Calcular masa adicional por nodos
nodes_dict = {}
for n, xyz in zip(nodes_tags, nodes_xyz):
    nodes_dict[n] = [xyz, 0.]   # coord, masa_adicional

thermal_elements = {}


tol = 1e-6
check_fixities = True
#Revisar que el modelo este apoyado en la caja de 6.6 m x 2.6 m x 7.8 m
for fixity in nodes_fixities:
    n = fixity[0]
    x = nodes_dict[n][0]
    # check_fixities &= abs(x[0] - Lx/2) <= tol
    # check_fixities &= abs(x[1] - Ly/2) <= tol
    # check_fixities &= abs(x[2] - Lz/2) <= tol

    if abs(x[0])-Lx/2 <= tol:
        if abs(x[1]) <= Ly/2 and abs(x[2]) <= Lz/2:
            pass
        else:
            check_fixities = False
            break
    if abs(x[1])-Ly/2 <= tol:
        if abs(x[0]) <= Lx/2 and abs(x[2]) <= Lz/2:
            pass
        else:
            check_fixities = False
            break
    if abs(x[2])-Lz/2 <= tol:
        if abs(x[0]) <= Lx/2 and abs(x[1]) <= Ly/2:
            pass
        else:
            check_fixities = False
            break

print(f"Cumple soportes       : {check_fixities}")
print("")


#Calcular la superficie y masa del panel
m_total_panel = 0.
inercia_total_panel = 0.
A_total_panel = 0.
for n in panel_nodes:
    x0 = nodes_dict[n[0]][0]
    x1 = nodes_dict[n[1]][0]
    x2 = nodes_dict[n[2]][0]

    d1 = x1 - x0
    d2 = x2 - x0
    d1xd2 = np.cross(d1,d2)/2

    A = np.sqrt(np.dot(d1xd2,d1xd2))
    vector_normal = d1xd2 / A

    m_panel_este_triangulo = A*ρ_panel

    A_total_panel += A
    m_total_panel += m_panel_este_triangulo
 
    nodes_dict[n[0]][1] += m_panel_este_triangulo/3
    nodes_dict[n[1]][1] += m_panel_este_triangulo/3
    nodes_dict[n[2]][1] += m_panel_este_triangulo/3
    
    inercia_total_panel += m_panel_este_triangulo/3 * (x0**2 + x1**2 + x2**2).sum()

    thermal_elements[(n[0], n[1])] = True
    thermal_elements[(n[0], n[2])] = True
    thermal_elements[(n[1], n[2])] = True

check_area = A_total_panel > area_min

print(f"Area total de panel   : {A_total_panel} m²")
print(f"Area requerida        : {area_min} m²")
print(f"Suficiente area       : {check_area}")
print("")

#Calcular la masa estructural
m_total_estructura = 0.
inercia_total_estructura = 0.
for e, nodos, props in zip(elements_tags, elements_connectivities, elements_section_info):
    Dint = props[0]
    t = props[1]
    A = np.pi * ((Dint/2 + t)**2 - (Dint/2)**2 )
    ρ_lin = A * ρ_compuesto
    x0 = nodes_dict[nodos[0]][0]
    x1 = nodes_dict[nodos[1]][0]
    L = np.sqrt(np.dot(x0-x1, x0-x1))
    m_total_estructura += L*ρ_lin

    inercia_total_estructura += L*ρ_lin * (((x0+x1)/2)**2).sum()


RME = m_total_estructura / m_total_panel * 100

check_RME = RME < RME_max

print(f"Masa total panel      : {m_total_panel} kg")
print(f"Masa total estructura : {m_total_estructura} kg")
print(f"Iner total panel      : {inercia_total_panel} kg*m^2")
print(f"Iner total estructura : {inercia_total_estructura} kg*m^2")
print(f"RME                   : {RME} %")
print(f"Cumple RME            : {check_area}")
print("")

def setup_model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)


    for n, xyz in zip(nodes_tags, nodes_xyz):
        m = nodes_dict[n][1]
        ops.node(int(n), xyz[0], xyz[1], xyz[2], "-mass", m, m, m)

    # Step 3: Define Material
    thermal_strain = alpha * ΔT

    mat_tag = 1
    ops.uniaxialMaterial('Elastic', mat_tag, E_compuesto)  # Elastic material with E=3000


    mat_tag_thermal = 2
    ops.uniaxialMaterial('InitStrainMaterial', mat_tag_thermal, mat_tag, -thermal_strain)  # material con deformacion inicial


    m_total_estructura = 0.

    for e, nodos, props in zip(elements_tags, elements_connectivities, elements_section_info):
        Dint = props[0]
        t = props[1]
        A = np.pi * ((Dint/2 + t)**2 - (Dint/2)**2 )
        ρ_lin = A * ρ_compuesto


        mtag = mat_tag
        if (nodos[0], nodos[1]) in thermal_elements \
            or (nodos[1], nodos[0]) in thermal_elements:
            mtag = mat_tag_thermal
        ops.element('Truss', int(e), int(nodos[0]), int(nodos[1]), A, mtag, "-rho", ρ_lin)

    # Apply boundary conditions (fix nodes at the base)
    for fixity in nodes_fixities:
        n = int(fixity[0])
        rx = int(fixity[1])
        ry = int(fixity[2])
        rz = int(fixity[3])
        ops.fix(n, rx, ry, rz)

def analyze_static(acc_x=0., acc_y=0., acc_z=0., ΔT=0.):

    setup_model()

    tsTag = 1  #identificador
    ops.timeSeries('Constant', tsTag)

    patternTag_x = 1
    patternTag_y = 2
    patternTag_z = 3
    dir_x = 1    # 1=x  2=y  3=z
    dir_y = 2    # 1=x  2=y  3=z
    dir_z = 3    # 1=x  2=y  3=z

    if acc_x != 0.:
        ops.pattern('UniformExcitation', patternTag_x, dir_x, "-accel", tsTag, "-fact", acc_x)
    if acc_y != 0.:
        ops.pattern('UniformExcitation', patternTag_y, dir_y, "-accel", tsTag, "-fact", acc_y)
    if acc_z != 0.:
        ops.pattern('UniformExcitation', patternTag_z, dir_z, "-accel", tsTag, "-fact", acc_z)


    factor_de_carga = 1.0
    ops.system('BandGeneral')        # sistema de ecuaciones 
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('LoadControl', factor_de_carga)   # analisis estatico carga controlada
    ops.algorithm('Linear')              # problema lineal  --> Newton pa no-lineal
    ops.analysis('Static')               # estatico 

    numero_de_pasos_de_analisis = 1
    ops.analyze(numero_de_pasos_de_analisis)

    desplazamientos = {n: ops.nodeDisp(n)[0] for n in ops.getNodeTags()}
    fuerzas_axiales = {e: ops.eleResponse(e, "axialForce")[0] for e in ops.getEleTags()}

    ops.wipe()

    return desplazamientos, fuerzas_axiales

def analyze_eigen():

    setup_model()

    # Step 5: Eigenvalue Analysis
    num_modes = 10
    eigenvalues = ops.eigen(num_modes)    #K phi = w² M phi

    # Calculate natural frequencies (Hz)
    eigenfrequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
    # print("Eigenfrequencies: ", eigenfrequencies)

    ops.wipe()

    return eigenfrequencies[0]





#Check fuerzas
disp_X, f_X = analyze_static(acc_x = 1)
disp_Y, f_Y = analyze_static(acc_y = 1)
disp_Z, f_Z = analyze_static(acc_y = 1)
disp_thermal, f_thermal = analyze_static( ΔT=ΔT)

F = {e: a_max * np.sqrt(f_X[e]**2 + f_Y[e]**2 + f_Z[e]**2) for e in f_X}

check_resistencia = True
check_pandeo = True
FU_min = np.inf
FU_max = -np.inf
for e, nodos, props in zip(elements_tags, elements_connectivities, elements_section_info):
    Dint = props[0]
    t = props[1]
    A = np.pi * ((Dint/2 + t)**2 - (Dint/2)**2 )
    I = np.pi * ((Dint/2 + t)**4 - (Dint/2)**4 )

    f = F[e]

    x0 = nodes_dict[nodos[0]][0]
    x1 = nodes_dict[nodos[1]][0]
    L = np.sqrt(np.dot(x0-x1, x0-x1))

    # check_pandeo
    Pcr = np.pi**2 * E_compuesto * I / (L**2)

    check_resistencia = True
    check_pandeo = True

    check_pandeo_this = FS*F[e] < Pcr
    if not check_pandeo_this:
        check_pandeo = False

    # check resistencia
    F_ultimate = σ_max_compuesto * A
    check_resistencia_this = FS*F[e] < F_ultimate

    if not check_resistencia_this:
        check_resistencia = False


    # print(f"{e = } {FS*F[e]= } {Pcr=} {F_ultimate=} {check_pandeo_this=} {check_resistencia_this=}")


    FU = FS * F[e] / F_ultimate

    FU_max = max(FU, FU_max)
    FU_min = min(FU, FU_min)


print(f"Cumple resistencia    : {check_resistencia}   {FU_min} < FU < {FU_max}")
print(f"Cumple pandeo         : {check_pandeo}")
print(f"")

f1 = analyze_eigen()
check_frequencia = f1 > f_min_requerida
print(f"Frecuencia fundamental: {f1} Hz")
print(f"Cumple f1             : {check_frequencia}")


θ_all_max = 0.
for n in panel_nodes:
    x0 = nodes_dict[n[0]][0]
    x1 = nodes_dict[n[1]][0]
    x2 = nodes_dict[n[2]][0]

    #Calcular la normal no-deformada
    d1 = x1 - x0
    d2 = x2 - x0
    d1xd2 = np.cross(d1,d2)
    n1 = d1xd2 / np.sqrt(np.dot(d1xd2,d1xd2))

    #Calcular la normal deformada
    x0 += disp_thermal[n[0]]
    x1 += disp_thermal[n[1]]
    x2 += disp_thermal[n[2]]
    d1 = x1 - x0
    d2 = x2 - x0
    d1xd2 = np.cross(d1,d2)
    n2 = d1xd2 / np.sqrt(np.dot(d1xd2,d1xd2))

    #rotacion
    cosθ = np.dot(n1, n2)
    θ = abs(np.rad2deg(np.arccos(cosθ)))
    θ_all_max = max(θ_all_max, θ)


check_θ_max = θ_all_max < θ_max

print(f"θ_all_max             : {θ_all_max} °")
print(f"Cumple θ_max          : {check_θ_max}")
print("")
# print(f"{f_X=}")
# print(f"{f_Y=}")
# print(f"{f_Z=}")
# print(f"{f_thermal=}")
# print(f"{disp_thermal=}")
# print(f"{f1=}")