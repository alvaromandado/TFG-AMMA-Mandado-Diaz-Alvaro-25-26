# GRAPE and CRAB Algorithms with SciPy

*Read this in other languages: [Español](README.md)*

This directory contains the manual implementation of the GRAPE algorithm (L-BFGS-B) and the CRAB algorithm (Nelder-Mead) using the SciPy package. These are used to optimally realize the CNOT gate using the following Hamiltonian:
$$H=J(\sigma_z\otimes\sigma_z)+u_0(t) \cdot 1+(u_1(t)\sigma_x+u_2(t)\sigma_y)\otimes 1+1\otimes(u_3(t)\sigma_x+u_4(t)\sigma_y)$$

## Description
The file `GRAPE_SciPy.py` defines the Hamiltonian and problem parameters, then optimizes the controls via *scipy.optimize.minimize*.
Optimization is performed by running the algorithm with random initial controls, evaluating the objective function, and displaying the controls that yield the best optimization.

The file `CRAB_SciPy.py` defines the Hamiltonian and problem parameters, then optimizes the controls via *scipy.optimize.minimize*.
Optimization is performed by running the algorithm for random frequencies, evaluating the parameters to be optimized, and displaying the controls that yield the best optimization.

## Contents
* **`GRAPE_SciPy.py`**: Implementation of GRAPE using SciPy. The code runs the algorithm multiple times and selects the best result.
* **`grape_scipy.png`**: Graphical representation of the best controls found by the code for GRAPE.
* **`CRAB_SciPy.py`**: Implementation of CRAB using SciPy. The code runs the algorithm multiple times and selects the best result.
* **`crab_scipy.png`**: Graphical representation of the best controls found by the code for CRAB.

## Execution
```bash
python GRAPE_SciPy.py
python CRAB_SciPy.py
