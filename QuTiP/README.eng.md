# GRAPE and CRAB Algorithms with QuTiP

*Read this in other languages: [Español](README.md)*

This directory contains the implementation of the GRAPE and CRAB algorithms using the QuTiP package. These are used to optimally realize the CNOT gate using the following Hamiltonian:
$$H=J(\sigma_z\otimes\sigma_z)+u_0(t) \cdot 1+(u_1(t)\sigma_x+u_2(t)\sigma_y)\otimes 1+1\otimes(u_3(t)\sigma_x+u_4(t)\sigma_y)$$

## Description
The file `GRAPE_QuTiP.py` describes the Hamiltonian and problem parameters to then optimize the controls through *qutip*.
Optimization is performed by running the algorithm with random initial controls, evaluating the objective function values, and displaying the controls that best optimize the function.

The file `CRAB_QuTiP.py` describes the Hamiltonian and problem parameters to then optimize the controls through *qutip*.
Optimization is performed by running the algorithm for random frequencies, evaluating the parameters to be optimized, and displaying the controls that best optimize the function.

## Contents
* **`GRAPE_QuTiP.py`**: Implementation of GRAPE using QuTiP. The code runs the algorithm multiple times and selects the best result.
* **`grape_qutip.png`**: Graphical representation of the best controls found by the code for GRAPE.
* **`CRAB_QuTiP.py`**: Implementation of CRAB using QuTiP. The code runs the algorithm multiple times and selects the best result.
* **`crab_qutip.png`**: Graphical representation of the best controls found by the code for CRAB.

## Execution
```bash
python GRAPE_QuTiP.py
python CRAB_QuTiP.py
