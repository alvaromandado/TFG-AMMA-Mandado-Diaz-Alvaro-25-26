# Algoritmos GRAPE y CRAB con QuTiP

Este directorio contiene la implementación manual del algoritmo GRAPE (L-BFGS-B) y del algoritmo CRAB (Nelder-Mead) a través del paquete de QuTiP. Este se utiliza para alcanzar óptimamente la puerta CNOT a partir del siguiente hamiltoniano:
$$H=J(\sigma_z\otimes\sigma_z)+u_0(t) \cdot 1+(u_1(t)\sigma_x+u_2(t)\sigma_y)\otimes 1+1\otimes(u_3(t)\sigma_x+u_4(t)\sigma_y)$$

## Descripción
El archivo `GRAPE_QuTiP.py` describe el hamiltoniano y los parámetros del problema para luego optimizar los controles a través de *qutip*.
Se optimiza corriendo el algoritmo con controles iniciales aleatorios, se calculan los valores de la función a optimizar y se muestran los controles que más optimizan la función.

El archivo `CRAB_QuTiP.py` describe el hamiltoniano y los parámetros del problema para luego optimizar los controles a través de *qutip*.
Se optimiza corriendo el algoritmo para frecuencias aleatorias, se calculan los valores de los parámetros a optimizar y se muestran los controles que más optimizan la función.

## Contenido
* **`GRAPE_SciPy.py`**: Implementación de GRAPE utilizando QuTiP. El código corre el algoritmo varias veces y elige el mejor resultado.
* **`grape_qutip.png`**: Representación gráfica de los mejores controles encontrados por el código para GRAPE.
* **`CRAB_QuTiP.py`**: Implementación de CRAB utilizando QuTiP. El código corre el algoritmo varias veces y elige el mejor resultado.
* **`crab_qutip.png`**: Representación gráfica de los mejores controles encontrados por el código para CRAB.


## Ejecución
```bash
python GRAPE_QuTiP.py
python CRAB_QuTiP.py
```
