# Algoritmos GRAPE y CRAB con SciPy

Este directorio contiene la implementación manual del algoritmo GRAPE (L-BFGS-B) y del algoritmo CRAB (Nelder-Mead) a través del paquete de SciPy. Este se utiliza para alcanzar óptimamente la puerta CNOT a partir del siguiente hamiltoniano:
$$H=J(\sigma_z\otimes\sigma_z)+(u_1(t)\sigma_x+u_2(t)\sigma_y)\otimes 1+1\otimes(u_3(t)\sigma_x+u_4(t)\sigma_y)$$

## Descripción
El archivo `GRAPE_SciPy.py` describe el hamiltoniano y los parámetros del problema para luego optimizar los controles a través de *scipy.optimize.minimize*.
Se optimiza corriendo el algoritmo con controles iniciales aleatorios, se calculan los valores de la función a optimizar y se muestran los controles que más mejoran la función.

El archivo `CRAB_SciPy.py` describe el hamiltoniano y los parámetros del problema para luego optimizar los controles a través de *scipy.optimize.minimize*.
Se optimiza corriendo el algoritmo con frecuencias aleatorias, se calculan los valores de los parámetros a optimizar y se muestran los controles que más mejoran la función.

## Contenido
* **`GRAPE_SciPy.py`**: Implementación de GRAPE utilizando SciPy. El código corre el algoritmo varias veces y elige el mejor resultado.
* **`grape_scipy.png`**: Representación gráfica de los mejores controles encontrados por el código para GRAPE.
* **`CRAB_SciPy.py`**: Implementación de CRAB utilizando SciPy. El código corre el algoritmo varias veces y elige el mejor resultado.
* **`crab_scipy.png`**: Representación gráfica de los mejores controles encontrados por el código para CRAB.

## Ejecución
```bash
python GRAPE_SciPy.py
python CRAB_SciPy.py
```
