# Appendix: Introduction to Optimal Control in the Dynamics of Quantum Systems (Department of Mathematical Analysis and Applied Mathematics)

*Read this in other languages: [Español](README.md)*

This repository contains the source code and generated figures for the Bachelor's Thesis titled **"Introduction to Optimal Control in the Dynamics of Quantum Systems"**, presented at the Faculty of Mathematical Sciences of the Complutense University of Madrid (Academic Year 2025-2026).

## Files and Results Correspondence

The following table links the sections of the memoir with their technical implementation and the graphical results obtained:

| Section | Directory | Code File | Generated Figure | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Chap. 3** | `/QuTiP` | `GRAPE_QuTiP.py`, `CRAB_QuTiP.py` | `grape_qutip.png`, `crab_qutip.png` | Use of the QuTiP package for optimal CNOT gate search. |
| **Chap. 3** | `/SciPy` | `GRAPE_SciPy.py`, `CRAB_SciPy.py` | `grape_scipy.png`, `crab_scipy.png` | Use of the SciPy package for optimal CNOT gate search. |

## Prerequisites and Installation

The simulations have been developed in **Python 3.10+**. To run them, it is recommended to create a virtual environment and install the necessary dependencies:

```bash
pip install numpy matplotlib qutip qutip-qtrl qutip-qoc joblib scipy
