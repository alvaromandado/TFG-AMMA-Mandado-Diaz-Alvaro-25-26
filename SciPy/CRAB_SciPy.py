# CÓDIGO PARA EL CÁLCULO NUMÉRICO CON CRAB Y EL PAQUETE SCIPY
# Este código forma parte de un anexo adicional del Trabajo de Fin de Grado para la obtención del Doble Grado en Física y Matemáticas de nombre
# "Introducción a la teoría de control óptimo en la dinámica de sistemas cuánticos", del curso 2026. Los tutores son Mihaela Negreanu Pruna y Federico Herrero Hervás.
# El autor de este código, así como del Trabajo de Fin de Grado al que pertenece, es Álvaro Mandado Díaz.

# Este código incluye la implementación del algoritmo CRAB con QuTiP para resolver el problema de la puerta CNOT. Se irá comentando qué hace cada parte del mismo.

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
dt = t_f / N

# Parámetros específicos de CRAB
n_controles = 5
n_frecuencias = 3   # Número de armónicos por control (3 senos y 3 cosenos)
# Total de variables: 5 controles * (3 A + 3 B + 1 C) = 35 variables

# Tolerancias y límites de Nelder-Mead
tol_error = 1e-4
max_iteraciones = 3000

# Definimos las puertas lógicas a utilizar
I = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

I_4 = np.eye(4, dtype=complex)

# Hamiltoniano de deriva original
H_d = J_const * np.kron(sz, sz)

# Hamiltonianos de control (u0, u1, u2, u3, u4)
H_c = [
    I_4,
    np.kron(sx, I), np.kron(sy, I),
    np.kron(I, sx), np.kron(I, sy)
]

# Matriz objetivo (CNOT)
U_target = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
], dtype=complex)

tiempos = np.linspace(0, t_f, N)

# Definimos la envolvente global: sin^2(pi * t / t_f)
envolvente = np.sin(np.pi * tiempos / t_f)**2

def ejecutar_intento_crab(intento_id, semilla):
    inicio_tiempo = time.time()
    np.random.seed(semilla)
    
    # 1. Definir las frecuencias base
    omegas = np.zeros((n_controles, n_frecuencias))
    for k in range(n_controles):
        for m in range(n_frecuencias):
            freq_fundamental = (2 * np.pi * (m + 1)) / t_f
            ruido = np.random.uniform(-0.2, 0.2)
            omegas[k, m] = freq_fundamental * (1 + ruido)

    # Matrices de senos y cosenos
    matriz_senos = np.zeros((n_controles, n_frecuencias, N))
    matriz_cosenos = np.zeros((n_controles, n_frecuencias, N))
    
    for k in range(n_controles):
        for m in range(n_frecuencias):
            matriz_senos[k, m, :] = np.sin(omegas[k, m] * tiempos)
            matriz_cosenos[k, m, :] = np.cos(omegas[k, m] * tiempos)

    # 2. Definimos la función de coste
    def evaluar_funcional(coeficientes_flat):
        # Desempaquetamos los 35 coeficientes
        A = coeficientes_flat[0:15].reshape((n_controles, n_frecuencias))
        B = coeficientes_flat[15:30].reshape((n_controles, n_frecuencias))
        C = coeficientes_flat[30:35]
        
        # Reconstruimos los controles añadiendo la constante C[k] y multiplicando por la envolvente
        u = np.zeros((n_controles, N))
        for k in range(n_controles):
            base_con_offset = C[k] + np.sum(A[k, :, np.newaxis] * matriz_senos[k, :, :] + 
                                            B[k, :, np.newaxis] * matriz_cosenos[k, :, :], axis=0)
            u[k, :] = envolvente * base_con_offset
            
        # Evolución temporal
        U = np.eye(4, dtype=complex)
        for j in range(N):
            H_total = H_d + sum(u[k, j] * H_c[k] for k in range(n_controles))
            U = expm(-1j * H_total * dt) @ U
            
        # Fidelidad Original (Parte real de la traza)
        g = -np.real(np.trace(U_target.conj().T @ U))
        
        # Fluencia
        f = (alpha / 2.0) * np.sum(u**2) * dt
        
        return g + f

    # 3. Optimización
    # Inicializamos los coeficientes cerca de cero
    coeficientes_iniciales = np.random.randn(n_controles * n_frecuencias * 2 + n_controles) * 0.1
    
    # Optimizador Nelder-Mead
    resultado = minimize(
        fun=evaluar_funcional,
        x0=coeficientes_iniciales,
        method='Nelder-Mead',
        options={
            'maxiter': max_iteraciones,
            'maxfev': max_iteraciones * 2,
            'disp': False, 
            'xatol': tol_error,
            'fatol': tol_error
        }
    )
    
    fin_tiempo = time.time()
    
    # Reconstruimos óptimos finales
    A_opt = resultado.x[0:15].reshape((n_controles, n_frecuencias))
    B_opt = resultado.x[15:30].reshape((n_controles, n_frecuencias))
    C_opt = resultado.x[30:35]
    
    u_opt = np.zeros((n_controles, N))
    for k in range(n_controles):
        base_con_offset = C_opt[k] + np.sum(A_opt[k, :, np.newaxis] * matriz_senos[k, :, :] + 
                                            B_opt[k, :, np.newaxis] * matriz_cosenos[k, :, :], axis=0)
        u_opt[k, :] = envolvente * base_con_offset
                             
    # Matriz final
    U_final_optimo = np.eye(4, dtype=complex)
    for j in range(N):
        H_total = H_d + sum(u_opt[k, j] * H_c[k] for k in range(n_controles))
        U_final_optimo = expm(-1j * H_total * dt) @ U_final_optimo

    return {
        'Intento': intento_id + 1,
        'Funcional_J': resultado.fun,
        'Tiempo_s': fin_tiempo - inicio_tiempo,
        'u_opt': u_opt,
        'U_final': U_final_optimo
    }

# Ejecutamos los intentos en paralelo
numero_intentos = 5

resultados_brutos = Parallel(n_jobs=-1)(
    delayed(ejecutar_intento_crab)(i, i+400) for i in range(numero_intentos)
)

# Procesamos los resultados para obtener el mejor resultado global
resultados_tabla = []
mejor_J = float('inf')
mejor_resultado = None

for res in resultados_brutos:
    resultados_tabla.append(res)
    if res['Funcional_J'] < mejor_J:
        mejor_J = res['Funcional_J']
        mejor_resultado = res

print("\n" + "="*65)
print("Resultados CRAB (Nelder-Mead)")
print("-" * 65)
print(f"{'Intento':<10} | {'Funcional J':<15} | {'Tiempo (s)':<10}")
print("-" * 65)
for fila in resultados_tabla:
    print(f"{fila['Intento']:<10} | {fila['Funcional_J']:<15.4f} | {fila['Tiempo_s']:<10.2f}")
print("="*65)


print("Matriz final obtenida:")
print(np.round(mejor_resultado['U_final'], decimals=3))
print("\nMatriz objetivo (teórica CNOT):")
print(np.round(U_target, decimals=3))

# --- GRÁFICAS ---
u_opt = mejor_resultado['u_opt']

fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

# Subplot 0: u0
max_u0 = np.max(np.abs(u_opt[0, :])) if np.max(np.abs(u_opt[0, :])) > 0 else 1.0
ax[0].plot(tiempos, envolvente * max_u0, '--', color='gray', alpha=0.4, label='Envolvente')
ax[0].plot(tiempos, -envolvente * max_u0, '--', color='gray', alpha=0.4)
ax[0].plot(tiempos, u_opt[0, :], label='$u_0(t)$, $\mathbb{I}\otimes\mathbb{I}$', color='#7f7f7f')
ax[0].set_ylabel('Amplitud')
ax[0].legend(loc='upper right')
ax[0].grid(True, alpha=0.3)
ax[0].set_title(f'Controles óptimos numéricos (CRAB con SciPy, intento: {mejor_resultado["Intento"]}) ($J$ = {mejor_J:.3f})')

# Subplot 1: u1 y u2
max_u12 = np.max(np.abs(u_opt[1:3, :]))
ax[1].plot(tiempos, envolvente * max_u12, '--', color='gray', alpha=0.4, label='Envolvente')
ax[1].plot(tiempos, -envolvente * max_u12, '--', color='gray', alpha=0.4)
ax[1].plot(tiempos, u_opt[1, :], label='$u_1(t)$, $\sigma_x\otimes1$', color='#1f77b4')
ax[1].plot(tiempos, u_opt[2, :], label='$u_2(t)$, $\sigma_y\otimes1$', color='#ff7f0e')
ax[1].set_ylabel('Amplitud')
ax[1].legend(loc='upper right')
ax[1].grid(True, alpha=0.3)

# Subplot 2: u3 y u4
max_u34 = np.max(np.abs(u_opt[3:5, :]))
ax[2].plot(tiempos, envolvente * max_u34, '--', color='gray', alpha=0.4, label='Envolvente')
ax[2].plot(tiempos, -envolvente * max_u34, '--', color='gray', alpha=0.4)
ax[2].plot(tiempos, u_opt[3, :], label='$u_3(t)$, $1\otimes\sigma_x$', color='#2ca02c')
ax[2].plot(tiempos, u_opt[4, :], label='$u_4(t)$, $1\otimes\sigma_y$', color='#d62728')
ax[2].set_xlabel('Tiempo ($t$)')
ax[2].set_ylabel('Amplitud')
ax[2].legend(loc='upper right')
ax[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
