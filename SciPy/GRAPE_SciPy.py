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

# Hamiltoniano de deriva, como operador producto

H_d_original = J_const * np.kron(sz, sz)

# Hamiltoniano de deriva actualizado con la fase particular

fase_adhoc = -np.pi / 4
H_d = H_d_original + (fase_adhoc / t_f) * I_4

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

# Definimos una función que calcula en cada iteración el gradiente y el valor de J, que es el rendimiento

def evaluar_funcional_y_gradiente(u_flat):
    
    # u_flat es un array 1D que L-BFGS-B maneja, lo reformateamos a  tamaño (4, N)
    
    u = u_flat.reshape((4, N))
    
    # El operador U evoluciona hacia adelante
    
    U = np.zeros((N + 1, 4, 4), dtype=complex)   # U(t_n)=U[n,:,:], en la primera coordenada se guarda el tiempo
    U[0] = np.eye(4, dtype=complex)
    
    for j in range(N):
        H_total = H_d + sum(u[k, j] * H_c[k] for k in range(4))
        U[j+1] = expm(-1j * H_total * dt) @ U[j]
        
    # Calculamos el valor del funcional J
    
    U_final = U[-1]
    g = -np.real(np.trace(U_target.conj().T @ U_final))
    f = (alpha / 2.0) * np.sum(u**2) * dt
    J_total = g + f
    
    # El coestado P evoluciona hacia atrás
    
    P = np.zeros((N + 1, 4, 4), dtype=complex)   # P(t_n)=P[n,:,:], en la primera coordenada se guarda el tiempo
    P[-1] = U_target # Condición de contorno del PMP
    
    for j in range(N - 1, -1, -1):   # Hacia atrás
        H_total = H_d + sum(u[k, j] * H_c[k] for k in range(4))
        P[j] = expm(1j * H_total * dt) @ P[j+1]   # El coestado evoluciona hacia atrás  (exponente positivo)
        
    # Calculamos el gradiente
    
    grad = np.zeros((4, N))
    for j in range(N):
        for k in range(4):
            # Derivada del término de superposición (usando aproximación de primer orden de expm)
            grad_g = -np.real(np.trace(P[j+1].conj().T @ (-1j * H_c[k] * dt) @ U[j]))   # Derivada de g
            grad_f = alpha * u[k, j] * dt   # Derivada explícita de f
            
            grad[k, j] = grad_g + grad_f
            
    return J_total, grad.flatten()   # Se aplana el vector para que las coordenadas i, i+1, i+2 e i+3 sean los controles u1, u2, u3 y u4 en la iteración i


# Sacamos los resultados por una tabla

resultados_tabla = []
mejor_J = float('inf')   # En minimización, el peor resultado posible es J = inf 
mejor_resultado = None
mejor_nombre_global = None

# Definimos las amplitudes máximas para todas las variables

limites = [(-5.0, 5.0)] * (4 * N)

# Hacemos 5 intentos con distintas semillas

intentos = 5
for i in range(intentos):

    inicio_tiempo = time.time()   # Para medir el tiempo de cálculo
    
    # Semilla inicial aleatoria en forma de ruido blanco pequeño
    u_inicial = np.random.randn(4 * N) * 0.1
    
    resultado = minimize(
        fun=evaluar_funcional_y_gradiente,   # Incluye los hamiltonianos y el cálculo de J
        x0=u_inicial,
        method='L-BFGS-B',   # Elegimos el método de optimización
        jac=True,
        bounds=limites,
        options={'maxiter': max_iteraciones, 'ftol': tol_error}
    )

    fin_tiempo = time.time()

    tiempo_ejecucion = fin_tiempo - inicio_tiempo

    resultados_tabla.append({
        'Intento': i+1,
        'Funcional_J': resultado.fun,
        'Tiempo_s': tiempo_ejecucion
    })

    # Guardamos el mejor intento de minimizar J para este generador
    
    if resultado.fun < mejor_J:
        mejor_J = resultado.fun
        mejor_resultado = resultado
        mejor_nombre_global = i+1

# Sacamos por la consola una tabla con los resultados

print("\n" + "="*85)
print(f"{'Intento':<15} | {'Funcional J':<15} | {'Tiempo (s)':<10}")
print("-" * 85)
for fila in resultados_tabla:
    print(f"{fila['Intento']:<15} | {fila['Funcional_J']:<15.4f} | {fila['Tiempo_s']:<10.2f}")
print("="*85)

# Sacamos por la consola el mejor resultado

print(f"\n=> El mejor conjunto de controles fue encontrado en el intento número {mejor_nombre_global} (J = {mejor_J:.4f})")


# Guardamos el mejor resultado

resultado = mejor_resultado

# Dibujamos los controles obtenidos

u_opt = resultado.x.reshape((4, N))
tiempos = np.linspace(0, t_f, N + 1)

fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

# Para el primer qubit dibujamos u_1 y u_2

ax[0].stairs(u_opt[0, :], tiempos, label='$u_1(t)$, $\sigma_x\otimes1$', color='#1f77b4', baseline=None)
ax[0].stairs(u_opt[1, :], tiempos, label='$u_2(t)$, $\sigma_y\otimes1$', color='#ff7f0e', baseline=None)
ax[0].set_ylabel('Amplitud')
ax[0].legend(loc='upper right')
ax[0].grid(True, alpha=0.3)
ax[0].set_title(f'Controles Óptimos (Funcional $J$ = {resultado.fun:.4f})')

# Para el segundo qubit dibujamos u_3 y u_4

ax[1].stairs(u_opt[2, :], tiempos, label='$u_3(t)$, $1\otimes\sigma_x$', color='#2ca02c', baseline=None)
ax[1].stairs(u_opt[3, :], tiempos, label='$u_4(t)$, $1\otimes\sigma_y$', color='#d62728', baseline=None)
ax[1].set_xlabel('Tiempo')
ax[1].set_ylabel('Amplitud')
ax[1].legend(loc='upper right')
ax[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# También podemos sacar por pantalla el operador obtenido

# Reconstruimos la evolución unitaria con los controles ganadores (u_opt)

U_final = np.eye(4, dtype=complex)
for j in range(N):
    H_total = H_d + sum(u_opt[k, j] * H_c[k] for k in range(4))
    U_final = expm(-1j * H_total * dt) @ U_final

# Redondeamos a tres decimales y sacamos por pantalla

matriz_limpia = np.round(U_final, decimals=3)
print(matriz_limpia)

# Comparamos con la matriz objetivo teórica

print(np.round(U_target, decimals=3))