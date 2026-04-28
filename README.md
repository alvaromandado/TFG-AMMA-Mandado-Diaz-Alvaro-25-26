# Anexo: Introducción al control óptimo en la dinámica de sistemas cuánticos (Departamento de Análisis Matemático y Matemática Aplicada)

Este repositorio contiene el código fuente y las figuras generadas para el Trabajo de Fin de Grado titulado **"Introducción al control óptimo en la dinámica de sistemas cuánticos"**, presentado en la Facultad de Ciencias Matemáticas de la Universidad Complutense de Madrid (Curso 2025-2026).

## Correspondencia de Archivos y Resultados

La siguiente tabla vincula las secciones de la memoria con su implementación técnica y los resultados gráficos obtenidos:

| Sección | Directorio | Archivo de Código | Figura Generada | Descripción |
| :--- | :--- | :--- | :--- | :--- |
| **Cap. 5.3**, Fig. 10 | `/QuTiP` | `GRAPE_QuTiP.py` | `qutip.png` | Uso del paquete QuTiP para encontrar óptimamente la puerta CNOT. |
| **Cap. 5.3**, Fig. 11 | `/SciPy` | `GRAPE_SciPy.py` | `scipy.png` | Uso de SciPy para escribir un algoritmo para encontrar óptimamente la puerta CNOT. |

## Requisitos e Instalación

Las simulaciones han sido desarrolladas en **Python 3.10+**. Para ejecutarlas, se recomienda crear un entorno virtual e instalar las dependencias necesarias:

```bash
pip install numpy matplotlib qutip qutip-qtrl qutip-qoc
```

### Instrucciones de uso:
1. Clone el repositorio: `git clone https://github.com/alvaromandado/TFG-AMMA-Mandado-Diaz-Alvaro-25-26.git`
2. Acceda a la carpeta del paquete deseado (ej: `cd QuTiP`).
3. Ejecute el script: `python GRAPE_QuTiP.py`.

## 📄 Licencia

Este proyecto está bajo la **Licencia MIT**. Siéntase libre de usar, modificar y distribuir el código, siempre que se mantenga la atribución al autor original.
