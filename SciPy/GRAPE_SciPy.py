# CÓDIGO PARA EL CÁLCULO NUMÉRICO CON GRAPE Y EL PAQUETE SCIPY
# Este código forma parte de un anexo adicional del Trabajo de Fin de Grado para la obtención del Doble Grado en Física y Matemáticas de nombre
# "Introducción a la teoría de control óptimo en la dinámica de sistemas cuánticos", del curso 2026. Los tutores son Mihaela Negreanu Pruna y Federico Herrero Hervás.
# El autor de este código, así como del Trabajo de Fin de Grado al que pertenece, es Álvaro Mandado Díaz.

# Este código incluye la implementación del algoritmo GRAPE con SciPy para resolver el problema de la puerta CNOT. Se irá comentando qué hace cada parte del mismo.

# Lo primero es importar los paquetes a utilizar.

import numpy as np    # Paquete para operaciones matemáticas
from scipy.linalg import expm   # Dentro de SciPy, sacar la exponenciación matricial
from scipy.optimize import minimize   # Dentro de SciPy, sacar la función de 'optimización'
import matplotlib.pyplot as plt    # Paquete para dibujar gráficas
import time    # Paquete para medir tiempo
from joblib import Parallel, delayed  # Paquete para paralelizar la ejecución de intentos de optimización

# Fijamos los parámetros del problema

J_const = 1.0       # Constante de acoplamiento ZZ
alpha = 0.05       # Penalización del uso de controles para el funcional J
t_f = 5.0           # Tiempo final
N = 1000            # Número de intervalos de la discretización
dt = t_f / N        # Paso temporal de la discretización
amp_limite = 5.0    # Límite de amplitud de los controles (para el algoritmo L-BFGS-B)

# Tolerancias y límites de la optimización

tol_error = 1e-4    # Tolerancia sobre el rendimiento
max_iteraciones = 10000    # Límite de seguridad de iteraciones
min_gradiente = 1e-8   # Tolerancia sobre el gradiente 

# Definimos las puertas lógicas a utilizar

I = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

I_4 = np.eye(4, dtype=complex)

# Hamiltoniano de deriva, como operador producto de dos qubits
H_d_original = J_const * np.kron(sz, sz)

# Añadimos los términos de control
H_c = [
    I_4,  # Hamiltoniano de control para u0 (Término de identidad global)
    np.kron(sx, I), np.kron(sy, I), # Hamiltonianos para u1 y u2
    np.kron(I, sx), np.kron(I, sy)  # Hamiltonianos para u3 y u4
]

# Definimos el operador evolución objetivo, la puerta CNOT
U_target = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
], dtype=complex)


# Puesto que L-BFGS-B requiere una función que devuelva el funcional y el gradiente a la vez,
# definimos la función principal que ejecutará la optimización GRAPE.

def optimizar_grape(H_drift, titulo_grafica):
    """
    Ejecuta el algoritmo GRAPE completo para un Hamiltoniano de deriva dado.
    """
    
    # Reiniciamos la semilla de números aleatorios para asegurar reproducibilidad.
    np.random.seed(0)

    # Definimos una función interna que calcula en cada iteración el gradiente y el valor de J
    def evaluar_funcional_y_gradiente(u_flat):
        
        # u_flat es un array 1D que L-BFGS-B maneja, lo reformateamos a tamaño (5, N)
        u = u_flat.reshape((5, N))
        
        # El operador U evoluciona hacia adelante
        U = np.zeros((N + 1, 4, 4), dtype=complex)   # U(t_n)=U[n,:,:], en la primera coordenada se guarda el tiempo
        U[0] = np.eye(4, dtype=complex)
        
        for j in range(N):
            H_total = H_drift + sum(u[k, j] * H_c[k] for k in range(5))
            U[j+1] = expm(-1j * H_total * dt) @ U[j]
            
        # Calculamos el valor del funcional J = g + f
        U_final = U[-1]
        g = -np.real(np.trace(U_target.conj().T @ U_final))
        f = (alpha / 2.0) * np.sum(u**2) * dt
        J_total = g + f
        
        # El coestado P evoluciona hacia atrás (Condición de contorno del PMP)
        P = np.zeros((N + 1, 4, 4), dtype=complex)   # P(t_n)=P[n,:,:]
        P[-1] = U_target 
        
        for j in range(N - 1, -1, -1):   # Propagación hacia atrás
            H_total = H_drift + sum(u[k, j] * H_c[k] for k in range(5))
            P[j] = expm(1j * H_total * dt) @ P[j+1]   # Evolución con exponente positivo
            
        # Calculamos el gradiente
        grad = np.zeros((5, N))
        for j in range(N):
            for k in range(5):
                # Derivada del término de superposición (usando aproximación de primer orden)
                grad_g = -np.real(np.trace(P[j+1].conj().T @ (-1j * H_c[k] * dt) @ U[j]))
                # Derivada explícita de la penalización
                grad_f = alpha * u[k, j] * dt
                
                grad[k, j] = grad_g + grad_f
                
        return J_total, grad.flatten()   # L-BFGS-B exige aplanar el vector del gradiente

    # --- INICIO DE LA PARALELIZACIÓN ---
    intentos = 5
    # Definimos las amplitudes máximas para todas las variables (ahora son 5 * N)
    limites = [(-amp_limite, amp_limite)] * (5 * N)

    # 1. Generamos los puntos de partida
    # Semilla inicial aleatoria en forma de ruido blanco pequeño
    u_iniciales = [np.random.randn(5 * N) * 0.1 for _ in range(intentos)]

    # Envolvemos un único intento en una función
    def ejecutar_intento(i):
        inicio_tiempo = time.time()   # Para medir el tiempo de cálculo
        
        u_inicial = u_iniciales[i]
        
        # Ejecución del algoritmo L-BFGS-B
        resultado = minimize(
            fun=evaluar_funcional_y_gradiente,   # Incluye los hamiltonianos y el cálculo de J
            x0=u_inicial,
            method='L-BFGS-B',   # Elegimos el método de optimización determinista
            jac=True,
            bounds=limites,
            options={'maxiter': max_iteraciones, 'ftol': tol_error, 'gtol': min_gradiente}
        )

        fin_tiempo = time.time()
        
        return {
            'Intento': i + 1,
            'Funcional_J': resultado.fun,
            'Tiempo_s': fin_tiempo - inicio_tiempo,
            'resultado_optimo': resultado
        }

    # Ahora ejecutamos los intentos en paralelo
    resultados_brutos = Parallel(n_jobs=-1)(delayed(ejecutar_intento)(i) for i in range(intentos))

    # Ahora procesamos los resultados para obtener el mejor resultado global
    resultados_tabla = []
    mejor_J = float('inf')   # En minimización, el peor resultado posible es J = inf 
    mejor_resultado = None
    mejor_nombre_global = None

    # Reconstruimos la tabla y buscamos el mejor resultado secuencialmente
    for res in resultados_brutos:
        resultados_tabla.append({
            'Intento': res['Intento'],
            'Funcional_J': res['Funcional_J'],
            'Tiempo_s': res['Tiempo_s']
        })
        
        # Guardamos el mejor intento de minimizar J
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
    u_opt = mejor_resultado.x.reshape((5, N))
    
    # Reconstruimos la evolución unitaria final con los controles ganadores (u_opt)
    U_final_optimo = np.eye(4, dtype=complex)
    for j in range(N):
        H_total = H_drift + sum(u_opt[k, j] * H_c[k] for k in range(5))
        U_final_optimo = expm(-1j * H_total * dt) @ U_final_optimo

    # Redondeamos a tres decimales y sacamos por pantalla
    print(f"Matriz final obtenida ({titulo_grafica}):")
    print(np.round(U_final_optimo, decimals=3))

    # Graficamos ese mejor resultado
    tiempos = np.linspace(0, t_f, N + 1)[:-1] # Ajuste para que coincida la dimensión de N con los controles
    
    # Ampliamos a 3 filas de subplots para incluir u0 de forma limpia
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    # Subplot 0: Para el término de u_0
    ax[0].stairs(u_opt[0, :], np.append(tiempos, t_f), label=r'$u_0(t)$, $1\otimes1$', color='#7f7f7f', baseline=None)
    ax[0].set_ylabel('Amplitud')
    ax[0].legend(loc='upper right')
    ax[0].grid(True, alpha=0.3)
    # Título adaptado para ser idéntico al del script de QuTiP
    ax[0].set_title(f'Controles óptimos numéricos (GRAPE con SciPy, intento: {mejor_nombre_global}) ($J$ = {mejor_J:.3f})')

    # Subplot 1: Para el primer qubit dibujamos u_1 y u_2
    ax[1].stairs(u_opt[1, :], np.append(tiempos, t_f), label=r'$u_1(t)$, $\sigma_x\otimes1$', color='#1f77b4', baseline=None)
    ax[1].stairs(u_opt[2, :], np.append(tiempos, t_f), label=r'$u_2(t)$, $\sigma_y\otimes1$', color='#ff7f0e', baseline=None)
    ax[1].set_ylabel('Amplitud')
    ax[1].legend(loc='upper right')
    ax[1].grid(True, alpha=0.3)

    # Subplot 2: Para el segundo qubit dibujamos u_3 y u_4
    ax[2].stairs(u_opt[3, :], np.append(tiempos, t_f), label=r'$u_3(t)$, $1\otimes\sigma_x$', color='#2ca02c', baseline=None)
    ax[2].stairs(u_opt[4, :], np.append(tiempos, t_f), label=r'$u_4(t)$, $1\otimes\sigma_y$', color='#d62728', baseline=None)
    ax[2].set_xlabel('Tiempo ($t$)')
    ax[2].set_ylabel('Amplitud')
    ax[2].legend(loc='upper right')
    ax[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

# Bloque principal protegido para la ejecución segura en paralelo en Windows
if __name__ == '__main__':
    # Ejecución para el Hamiltoniano original
    optimizar_grape(H_d_original, "Controles Óptimos GRAPE")

    # Comparamos con la matriz objetivo teórica
    print("\nMatriz objetivo (teórica CNOT):")
    print(np.round(U_target, decimals=3))
