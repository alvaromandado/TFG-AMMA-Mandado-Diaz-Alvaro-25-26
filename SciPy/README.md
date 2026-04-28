# Algoritmo GRAPE con SciPy

Este directorio contiene la implementación manual del algoritmo GRAPE (L-BFGS-B) a través del paquete de SciPy. Este se utiliza para alcanzar óptimamente la puerta CNOT a partir del siguiente hamiltoniano:
$$H=J(\sigma_z\otimes\sigma_z)+(u_1(t)\sigma_x+u_2(t)\sigma_y)\otimes 1+1\otimes(u_3(t)\sigma_x+u_4(t)\sigma_y)$$

## Descripción
El archivo `GRAPE_SciPy.py` describe el hamiltoniano y los parámetros del problema para luego optimizar los controles a través de *scipy.optimize.minimize*.
Se optimiza corriendo el algoritmo con controles iniciales aleatorios, se calculan los valores de la función a optimizar y se muestran los controles que más mejoran la función.

## Contenido
* **`GRAPE_SciPy.py`**: Implementación utilizando SciPy. El código corre el algoritmo varias veces y elige el mejor resultado.
* **`scipy.png`**: Representación gráfica de los mejores controles encontrados por el código.

## Ejecución
```bash
python GRAPE_SciPy.py
```
