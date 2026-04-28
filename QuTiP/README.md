# Algoritmo GRAPE con QuTiP

Este directorio contiene la implementación del algoritmo GRAPE (L-BFGS-B) a través del paquete de QuTiP. Este se utiliza para alcanzar óptimamente la puerta CNOT a partir del siguiente hamiltoniano:
$$H=J(\sigma_z\otimes\sigma_z)+(u_1(t)\sigma_x+u_2(t)\sigma_y)\otimes 1+1\otimes(u_3(t)\sigma_x+u_4(t)\sigma_y)$$

## Descripción
El archivo `GRAPE_QuTiP.py` describe el hamiltoniano y los parámetros del problema para luego optimizar los controles a través de *qutip.control.pulseoptim.optimize_pulse_unitary*.
Se optimiza mediante distintos generadores de pulsos iniciales, se calculan los valores de la función a optimizar y se muestran los controles que más mejoran la función.

## Contenido
* **`GRAPE_QuTiP.py`**: Implementación utilizando QuTiP. El código corre el algoritmo varias veces y elige el mejor resultado.
* **`qutip.png`**: Representación gráfica de los mejores controles encontrados por el código.

## Ejecución
```bash
python GRAPE_QuTiP.py
```
