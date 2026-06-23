# CÓDIGO PARA EL CÁLCULO NUMÉRICO CON CRAB Y EL PAQUETE SCIPY
# Este código forma parte de un anexo adicional del Trabajo de Fin de Grado para la obtención del Doble Grado en Física y Matemáticas de nombre
# "Introducción a la teoría de control óptimo en la dinámica de sistemas cuánticos", del curso 2026. Los tutores son Mihaela Negreanu Pruna y Federico Herrero Hervás.
# El autor de este código, así como del Trabajo de Fin de Grado al que pertence, es Álvaro Mandado Díaz.

# Este código incluye la implementación del algoritmo CRAB con SciPy para resolver el problema de la puerta CNOT. Se irá comentando qué hace cada parte del mismo.

# Lo primero es importar los paquetes a utilizar.

import numpy as np    # Paquete para operaciones matemáticas
from scipy.linalg import expm   # Dentro de SciPy, sacar la exponenciación matricial
from scipy.optimize import minimize   # Dentro de SciPy, sacar la función de 'optimización'
import matplotlib.pyplot as plt    # Paquete para dibujar gráficas
import time    # Paquete para medir tiempo
from joblib import Parallel, delayed  # Paquete para paralelizar la ejecución de intentos de optimización


# Fijamos los parámetros del problema
J_const = 1.0       # Constante de acoplamiento ZZ
alpha = 0.05        # Penalización del uso de controles para el funcional J
t_f = 5.0           # Tiempo final
N = 1000            # Número de intervalos de la discretización
dt = t_f / N        # Paso temporal de la discretización

# Parámetros específicos del algoritmo CRAB
N_c = 5             # Número de componentes de frecuencia por cada control
r_factor = 0.5      # Factor de aleatorización para las frecuencias

# Tolerancias y límites de la optimización
tol_error = 1e-4    # Tolerancia sobre el rendimiento
max_iteraciones = 10000    # Límite de seguridad de iteraciones

# Definimos las puertas lógicas a utilizar
I = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

I_4 = np.eye(4, dtype=complex)

# Hamiltoniano de deriva, como operador producto
H_d_original = J_const * np.kron(sz, sz)

# Hamiltonianos de control, como operadores producto
H_c = [
    np.kron(sx, I), np.kron(sy, I),
    np.kron(I, sx), np.kron(I, sy)
]

# Definimos el operador evolución objetivo, la puerta CNOT
U_target = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
], dtype=complex)

# Mallado temporal y  envolvente para CRAB
tiempos_eval = np.linspace(0, t_f, N, endpoint=False)
g_shape = np.sin(np.pi * tiempos_eval / t_f)**2    # Envolvente

def optimizar_crab(H_drift, titulo_grafica):
    """
    Ejecuta el algoritmo CRAB completo para un Hamiltoniano de deriva dado.
    """
    # Reiniciamos la semilla de números aleatorios para asegurar reproducibilidad. Aquí podría modificarse para obtener diferentes resultados en cada ejecución si se desea.
    np.random.seed(0)
    # Precalculamos las frecuencias aleatorizadas
    omegas = np.zeros((4, N_c))
    for i in range(4):
        for k in range(1, N_c + 1):
            rk = np.random.uniform(-r_factor, r_factor)
            omegas[i, k-1] = (2 * np.pi * k / t_f) * (1 + rk)

    # Definimos una función que calcula en cada iteración el valor de J a partir de las amplitudes
    def evaluar_funcional_J(c_flat):
        
        # c_flat es un array 1D que Nelder-Mead maneja, lo reformateamos a tamaño (4, N_c, 2)
        c = c_flat.reshape((4, N_c, 2))
        u = np.zeros((4, N))
        
        for i in range(4):
            for k in range(N_c):
                u[i, :] += c[i, k, 0] * np.sin(omegas[i, k] * tiempos_eval) + \
                           c[i, k, 1] * np.cos(omegas[i, k] * tiempos_eval)
        u *= g_shape    # Forzamos que empiece y termine en 0
        
        # El operador U evoluciona hacia adelante
        U_final = np.eye(4, dtype=complex)
        
        for j in range(N):
            H_total = H_drift + sum(u[k, j] * H_c[k] for k in range(4))
            U_final = expm(-1j * H_total * dt) @ U_final
            
        # Calculamos el valor del funcional J
        g_term = -np.real(np.trace(U_target.conj().T @ U_final))
        f_term = (alpha / 2.0) * np.sum(u**2) * dt
        
        return g_term + f_term

    # --- INICIO DE LA PARALELIZACIÓN ---
    intentos = 5
    dimension_c = 4 * N_c * 2

    # 1. Generamos los puntos de partida de forma SECUENCIAL para imitar tu código original
    # Esto usa la misma "cinta" iniciada por el np.random.seed(0) global
    c_iniciales = [np.random.randn(dimension_c) * 0.5 for _ in range(intentos)]

    # Envolvemos un único intento en una función
    def ejecutar_intento(i):
        inicio_tiempo = time.time()
        
        # 2. En lugar de generar aleatoriedad aquí dentro, usamos el array precalculado
        c_inicial = c_iniciales[i]
        
        resultado = minimize(
            fun=evaluar_funcional_J,
            x0=c_inicial,
            method='Nelder-Mead',
            options={'maxiter': max_iteraciones, 'xatol': tol_error, 'fatol': tol_error}
        )

        fin_tiempo = time.time()
        
        return {
            'Intento': i + 1,
            'Funcional_J': resultado.fun,
            'Tiempo_s': fin_tiempo - inicio_tiempo,
            'resultado_optimo': resultado
        }

    # Ejecutamos en paralelo
    resultados_brutos = Parallel(n_jobs=-1)(delayed(ejecutar_intento)(i) for i in range(intentos))

    # --- PROCESAMIENTO DE RESULTADOS ---
    resultados_tabla = []
    mejor_J = float('inf') 
    mejor_resultado = None
    mejor_nombre_global = None

    # Reconstruimos la tabla y buscamos el mejor resultado secuencialmente
    for res in resultados_brutos:
        resultados_tabla.append({
            'Intento': res['Intento'],
            'Funcional_J': res['Funcional_J'],
            'Tiempo_s': res['Tiempo_s']
        })
        
        if res['Funcional_J'] < mejor_J:
            mejor_J = res['Funcional_J']
            mejor_resultado = res['resultado_optimo']
            mejor_nombre_global = res['Intento']

    # Sacamos por la consola una tabla con los resultados
    print("\n" + "="*85)
    print(f"Resultados para: {titulo_grafica}")
    print("-" * 85)
    print(f"{'Intento':<15} | {'Funcional J':<15} | {'Tiempo (s)':<10}")
    print("-" * 85)
    for fila in resultados_tabla:
        print(f"{fila['Intento']:<15} | {fila['Funcional_J']:<15.4f} | {fila['Tiempo_s']:<10.2f}")
    print("="*85)

    # Sacamos por la consola el mejor resultado
    print(f"\n=> El mejor conjunto de controles fue encontrado en el intento número {mejor_nombre_global} (J = {mejor_J:.4f})\n")

    # Reconstruimos los controles óptimos a partir de las amplitudes
    c_opt = mejor_resultado.x.reshape((4, N_c, 2))
    u_opt = np.zeros((4, N))
    for i in range(4):
        for k in range(N_c):
            u_opt[i, :] += c_opt[i, k, 0] * np.sin(omegas[i, k] * tiempos_eval) + \
                           c_opt[i, k, 1] * np.cos(omegas[i, k] * tiempos_eval)
    u_opt *= g_shape

    # También podemos sacar por pantalla el operador obtenido
    # Reconstruimos la evolución unitaria con los controles ganadores (u_opt)
    U_final_optimo = np.eye(4, dtype=complex)
    for j in range(N):
        H_total = H_drift + sum(u_opt[k, j] * H_c[k] for k in range(4))
        U_final_optimo = expm(-1j * H_total * dt) @ U_final_optimo

    # Redondeamos a tres decimales y sacamos por pantalla
    print(f"Matriz final obtenida ({titulo_grafica}):")
    print(np.round(U_final_optimo, decimals=3))

    # Dibujamos los controles obtenidos
    tiempos = np.linspace(0, t_f, N + 1)
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    # Para el primer qubit dibujamos u_1 y u_2
    ax[0].stairs(u_opt[0, :], tiempos, label='$u_1(t)$, $\sigma_x\otimes1$', color='#1f77b4', baseline=None)
    ax[0].stairs(u_opt[1, :], tiempos, label='$u_2(t)$, $\sigma_y\otimes1$', color='#ff7f0e', baseline=None)
    ax[0].set_ylabel('Amplitud')
    ax[0].legend(loc='upper right')
    ax[0].grid(True, alpha=0.3)
    ax[0].set_title(f'{titulo_grafica} (Funcional $J$ = {mejor_resultado.fun:.4f})')

    # Para el segundo qubit dibujamos u_3 y u_4
    ax[1].stairs(u_opt[2, :], tiempos, label='$u_3(t)$, $1\otimes\sigma_x$', color='#2ca02c', baseline=None)
    ax[1].stairs(u_opt[3, :], tiempos, label='$u_4(t)$, $1\otimes\sigma_y$', color='#d62728', baseline=None)
    ax[1].set_xlabel('Tiempo')
    ax[1].set_ylabel('Amplitud')
    ax[1].legend(loc='upper right')
    ax[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# Ejecución para el Hamiltoniano original
optimizar_crab(H_d_original, "Controles Óptimos - H Original")

# Repetimos el proceso para el hamiltoniano modificado
fase_adhoc = -np.pi / 4
H_d_mod = H_d_original + fase_adhoc * I_4
optimizar_crab(H_d_mod, "Controles Óptimos - H Modificado")

# Comparamos con la matriz objetivo teórica
print("\nMatriz objetivo (teórica CNOT):")
print(np.round(U_target, decimals=3))
