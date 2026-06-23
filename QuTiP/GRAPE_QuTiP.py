# CÓDIGO PARA EL CÁLCULO NUMÉRICO CON GRAPE Y EL PAQUETE QUTIP
# Este código forma parte de un anexo adicional del Trabajo de Fin de Grado para la obtención del Doble Grado en Física y Matemáticas de nombre
# "Introducción a la teoría de control óptimo en la dinámica de sistemas cuánticos", del curso 2026. Los tutores son Mihaela Negreanu Pruna y Federico Herrero Hervás.
# El autor de este código, así como del Trabajo de Fin de Grado al que pertence, es Álvaro Mandado Díaz.

# Este código incluye la implementación del algoritmo GRAPE con QuTiP para resolver el problema de la puerta CNOT. Se irá comentando qué hace cada parte del mismo.

# Lo primero es importar los paquetes a utilizar.

import numpy as np    # Paquete para operaciones matemáticas
import matplotlib.pyplot as plt    # Paquete para dibujar gráficas
import time    # Paquete para medir tiempo
import qutip as qt    # Paquete de computación cuántica general
import qutip.control.pulseoptim as cpo    # Dentro de QuTiP, del paquete de control obtener el paquete de optimizar pulsos
from qutip.qip.operations import cnot    # Sacar la puerta CNOT, que será el objetivo
from joblib import Parallel, delayed  # Paquete para paralelizar la ejecución de intentos de optimización

# Fijamos los parámetros del problema

J_const = 1.0       # Constante de acoplamiento ZZ
alpha = 0.05        # Penalización del uso de controles para el funcional J
t_f = 5.0           # Tiempo final
N = 1000            # Número de intervalos de la discretización
dt = t_f / N
amp_limite = 5.0    # Límite de amplitud de los controles (para el algoritmo L-BFGS-B)

# Tolerancias y límites de la optimización
tol_error = 1e-4    # Tolerancia sobre el rendimiento
max_iteraciones = 10000    # Límite de seguridad de iteraciones
min_gradiente = 1e-8   # Tolerancia sobre el gradiente 

# Definimos las puertas lógicas a utilizar

I, sx, sy, sz = qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()    # Puertas de un qubit

H_d = J_const * qt.tensor(sz, sz)   # Hamiltoniano de deriva, como operador producto de dos qubits

H_c = [
    qt.tensor(sx, I), qt.tensor(sy, I),
    qt.tensor(I, sx), qt.tensor(I, sy)
]   # Hamiltonianos de control, como operadores producto de dos qubits

# Definimos los operadores evolución inicial y objetivo

U_0 = qt.tensor(I, I)
U_target = cnot()

# Puesto que L-BFGS-B maximiza el rendimiento pero no minimiza J, calculamos J explícitamente
# 'resultado' mide el final de la simulación en su conjunto

def funcional_J(resultado, alpha_val, dt_val, target):
    U_final = resultado.evo_full_final   # Obtenemos el operador evolución obtenido
    # g = -Re[tr(U_target^\dagger * U(t_f))]
    g_term = -np.real((target.dag() * U_final).tr())   # Calculamos el término g(U(t_f))
    controles = resultado.final_amps   # Importamos los controles finales
    penalizacion = (alpha_val / 2.0) * np.sum(controles**2) * dt_val   # Calculamos el termino de Lagrange
    
    return g_term + penalizacion   # Calculamos J = g + f

# Definimos los distintos tipos de generadores, estocásticos y deterministas

generadores_estocasticos = ['RND', 'RNDFOURIER', 'RNDWAVES']
generadores_deterministas = ['ZERO', 'SINE', 'LIN', 'SQUARE']

# Evaluamos todos los generadores, juntándolos en la misma lista

todos_los_generadores = generadores_estocasticos + generadores_deterministas


# Definimos una función interna para ejecutar un único intento en paralelo
def ejecutar_intento_qutip(gen, intento_id, semilla_proceso):
    inicio_tiempo = time.time()   # Para medir el tiempo de cálculo
    
    # Fijamos la semilla aquí dentro para obligar al proceso hijo a ser reproducible
    np.random.seed(semilla_proceso)
    
    # Ejecución del algoritmo GRAPE (en realidad L-BFGS-B)
    res = cpo.optimize_pulse_unitary(
        H_d, H_c, U_0, U_target,   # Hamiltonianos y condiciones de contorno
        num_tslots=N,   # Número de intervalos
        evo_time=t_f,   # Tiempo final
        amp_lbound=-amp_limite,   # Límites de amplitud
        amp_ubound=amp_limite,
        init_pulse_type=gen,   # Estimación inicial de pulsos de control
        fid_err_targ=tol_error,   # Tolerancia sobre la fidelidad
        max_iter=max_iteraciones,   # Límite superior de iteraciones
        min_grad=min_gradiente,   # Tolerancia sobre el gradiente
        # phase_option='SU',   # Fijar la fase, de forma que el elemento (1,1) de U(t_f) sea real positivo
        gen_stats=False
    )
    
    fin_tiempo = time.time()
    tiempo_ejecucion = fin_tiempo - inicio_tiempo
    
    # Calculamos Phi_0 y J
    phi_0 = 1.0 - res.fid_err  # Fidelidad normalizada a la unidad
    J_val = funcional_J(res, alpha, dt, U_target)   # Funcional J = g + f
    
    # Extraemos solo la información numérica estrictamente necesaria para aligerar la memoria en joblib
    final_amps_copia = np.array(res.final_amps)
    time_copia = np.array(res.time)
    u_final_copia = np.array(res.evo_full_final.full()) # Extraemos la matriz de evolución final como array de NumPy
    
    return {
        'Generador': gen,
        'Fidelidad_Phi0': phi_0,
        'Funcional_J': J_val,
        'Tiempo_s': tiempo_ejecucion,
        'final_amps': final_amps_copia,
        'time': time_copia,
        'U_final': u_final_copia
    }


# --- INICIO DE LA PARALELIZACIÓN ---

# Creamos la lista con la combinación de todas las tareas independientes a lanzar
# Asignamos un índice de semilla único e invariable a cada tarea individual
tareas = []
contador_semilla = 0
for gen in todos_los_generadores:
    es_estocastico = gen in generadores_estocasticos
    intentos = 5 if es_estocastico else 1
    for i in range(intentos):
        tareas.append((gen, i, contador_semilla))
        contador_semilla += 1

# Ejecutamos de forma masiva y en paralelo haciendo uso de todos los núcleos del procesador
resultados_brutos = Parallel(n_jobs=-1)(
    delayed(ejecutar_intento_qutip)(gen, i, sem) for gen, i, sem in tareas
)

# --- PROCESAMIENTO SECUENCIAL DE RESULTADOS ---

resultados_tabla = []
mejor_J_global = float('inf')   # En minimización, el peor resultado posible es J = inf
mejor_resultado_global = None
mejor_nombre_global = ""

# Reconstruimos los óptimos locales reuniendo ordenadamente los generadores
for gen in todos_los_generadores:
    es_estocastico = gen in generadores_estocasticos
    
    mejor_J_local = float('inf')   # En minimización, el peor resultado posible es J = inf
    mejor_phi_local = 0.0   # En términos de fidelidad, el peor resultado posible es Phi_0=0
    tiempo_local = 0.0
    resultado_local_optimo = None
    
    # Filtramos las ejecuciones correspondientes a este generador en particular
    intentos_del_generador = [r for r in resultados_brutos if r['Generador'] == gen]
    
    for res in intentos_del_generador:
        # Guardamos el mejor intento de minimizar J para este generador
        if res['Funcional_J'] < mejor_J_local:
            mejor_J_local = res['Funcional_J']
            mejor_phi_local = res['Fidelidad_Phi0']
            tiempo_local = res['Tiempo_s']
            resultado_local_optimo = {
                'time': res['time'],
                'final_amps': res['final_amps'],
                'U_final': res['U_final']
            }
            
    # Guardamos los resultados en la tabla
    resultados_tabla.append({
        'Generador': gen,
        'Tipo': 'Estocástico' if es_estocastico else 'Determinista',
        'Fidelidad_Phi0': mejor_phi_local,
        'Funcional_J': mejor_J_local,
        'Tiempo_s': tiempo_local
    })
    
    # Actualizamos el óptimo global para después dibujarlo
    if mejor_J_local < mejor_J_global:
        mejor_J_global = mejor_J_local
        mejor_resultado_global = resultado_local_optimo
        mejor_nombre_global = gen

# Sacamos por la consola una tabla con los resultados

print("\n" + "="*85)
print(f"{'Generador':<15} | {'Tipo':<13} | {'Fidelidad (Φ_0)':<17} | {'Funcional J':<15} | {'Tiempo (s)':<10}")
print("-" * 85)
for fila in resultados_tabla:
    print(f"{fila['Generador']:<15} | {fila['Tipo']:<13} | {fila['Fidelidad_Phi0']:<17.6f} | {fila['Funcional_J']:<15.4f} | {fila['Tiempo_s']:<10.2f}")
print("="*85)

# Sacamos por la consola el mejor resultado

print(f"\n=> El mejor conjunto de controles fue encontrado por: {mejor_nombre_global} (J = {mejor_J_global:.4f})")

# También podemos sacar por pantalla el operador obtenido

# Redondeamos a tres decimales y sacamos por pantalla

matriz_limpia = np.round(mejor_resultado_global['U_final'], decimals=3)
print(matriz_limpia)

# Comparamos con la matriz objetivo teórica

print(np.round(U_target.full(), decimals=3))

# Graficamos ese mejor resultado

tiempos = mejor_resultado_global['time']
u_opt = mejor_resultado_global['final_amps']

fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

# Para el primer qubit dibujamos u_1 y u_2

ax[0].stairs(u_opt[:, 0], tiempos, label='$u_1(t)$, $\sigma_x\otimes1$', color='#1f77b4', baseline=None)
ax[0].stairs(u_opt[:, 1], tiempos, label='$u_2(t)$, $\sigma_y\otimes1$', color='#ff7f0e', baseline=None)
ax[0].set_ylabel('Amplitud')
ax[0].legend(loc='upper right')
ax[0].grid(True, alpha=0.3)
ax[0].set_title(f'Mejores Controles - Generador: {mejor_nombre_global} ($J$ = {mejor_J_global:.3f})')

# Para el segundo qubit dibujamos u_3 y u_4

ax[1].stairs(u_opt[:, 2], tiempos, label='$u_3(t)$, $1\otimes\sigma_x$', color='#2ca02c', baseline=None)
ax[1].stairs(u_opt[:, 3], tiempos, label='$u_4(t)$, $1\otimes\sigma_y$', color='#d62728', baseline=None)
ax[1].set_xlabel('Tiempo ($t$)')
ax[1].set_ylabel('Amplitud')
ax[1].legend(loc='upper right')
ax[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
