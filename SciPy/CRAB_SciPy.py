import numpy as np    # Paquete para operaciones matemáticas
from scipy.linalg import expm   # Dentro de SciPy, sacar la exponenciación matricial
from scipy.optimize import minimize   # Dentro de SciPy, sacar la función de 'optimización'
import matplotlib.pyplot as plt    # Paquete para dibujar gráficas
import time    # Paquete para medir tiempo

# Para los casos estocásticos, fijamos una semilla para los paquetes de números pseudoaleatorios que permita la reproducibilidad 
np.random.seed(0)

# Fijamos los parámetros del problema
J_const = 1.0       # Constante de acoplamiento ZZ
alpha = 0.05        # Penalización del uso de controles para el funcional J
t_f = 5.0           # Tiempo final
N = 1000            # Número de intervalos de la discretización
dt = t_f / N

# Parámetros específicos del algoritmo CRAB
N_c = 5             # Número de componentes de frecuencia por cada campo de control
r_factor = 0.5      # Factor de aleatorización máxima de las frecuencias

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

# Tiempos y función de forma de CRAB
tiempos_eval = np.linspace(0, t_f, N, endpoint=False)
g_shape = np.sin(np.pi * tiempos_eval / t_f)**2    # Función para forzar contornos suaves

def optimizar_crab(H_drift, titulo_grafica):
    """
    Ejecuta el algoritmo CRAB completo para un Hamiltoniano de deriva dado.
    """
    
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

    # Sacamos los resultados por una tabla
    resultados_tabla = []
    mejor_J = float('inf')   # En minimización, el peor resultado posible es J = inf 
    mejor_resultado = None
    mejor_nombre_global = None

    # Hacemos 5 intentos con distintas semillas de amplitudes
    # Este proceso podría realizarse en paralelo en varios hilos mediante la librería joblib
    
    intentos = 5
    dimension_c = 4 * N_c * 2

    for i in range(intentos):
        inicio_tiempo = time.time()   # Para medir el tiempo de cálculo
        
        # Semilla inicial aleatoria para los coeficientes trigonométricos
        c_inicial = np.random.randn(dimension_c) * 0.5
        
        resultado = minimize(
            fun=evaluar_funcional_J,   # Solo evalúa J, no devuelve gradiente
            x0=c_inicial,
            method='Nelder-Mead',   # Elegimos el método del símplex clásico
            options={'maxiter': max_iteraciones, 'xatol': tol_error, 'fatol': tol_error}
        )

        fin_tiempo = time.time()
        tiempo_ejecucion = fin_tiempo - inicio_tiempo

        resultados_tabla.append({
            'Intento': i+1,
            'Funcional_J': resultado.fun,
            'Tiempo_s': tiempo_ejecucion
        })
        
        # Guardamos el mejor intento de minimizar J
        if resultado.fun < mejor_J:
            mejor_J = resultado.fun
            mejor_resultado = resultado
            mejor_nombre_global = i+1

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

# -----------------------------------------------------------------------------------
# EJECUCIÓN DEL CÓDIGO
# -----------------------------------------------------------------------------------

# 1. Ejecución para el Hamiltoniano libre original
optimizar_crab(H_d_original, "Controles Óptimos - H Libre")

# 2. Repetimos el proceso para el hamiltoniano modificado
fase_adhoc = -np.pi / 4
H_d_mod = H_d_original + fase_adhoc * I_4
optimizar_crab(H_d_mod, "Controles Óptimos - H Modificado")

# Comparamos con la matriz objetivo teórica
print("\nMatriz objetivo (teórica CNOT):")
print(np.round(U_target, decimals=3))
