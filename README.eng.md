# Appendix: Quantum Information and Quantum Computation (FT37)

*Read this in other languages: [Español](README.md)*


This repository contains the codes and figures released for the Bachelor's thesis **Quantum Information and Quantum Computation**, presented in Facultad de Ciencias Físicas at Universidad Complutense de Madrid (2025-2026)


## Files and results correspondance

The next table shows how the different sections in the memoir relate to the technical implementation and the graphic results:

| Section | Branch | Code file | Generated figure | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Sec. 2.2**, Fig. 2 | `/Deutsch-Jozsa` | `deutsch-jozsa-eng.py` | `deutsch-jozsa-eng.png` | Deutsch-Jozsa algorithm simulation, compared to the postquantum non-deterministic algorithm and random search. |
| **Sec. 2.4**, Fig. 3 | `/Grover` | `grover-eng.py` | `grover-eng.png` | Amplitud amplifying and randomized database search compaering Gorver's algorithm and linear search. |
| **Sec. 2.5**, Fig. 4 | `/Shor` | `shor-eng.py` | `shor-eng.png` | Integer factorization comparing Shor's algorithm with a series of classical algorithms. |

## Prerrequisites and Installation

The simulations have been implemented using **Python 3.10+**. To run them, it is highly recommended to use a virtual environment and to install the necessary packages:

```bash
pip install numpy matplotlib qiskit qiskit-aer
```

### Usage instructions:
1. Clone the repository: `git clone https://github.com/alvaromandado/TFG-FT37-Mandado-Diaz-Alvaro-25-26.git`
2. Access the desired algorithm's branch (e.g.: `cd Grover`).
3. Run the script: `python grover.py`.

## 📄 License

This project runs under a **MIT License**. Feel free to use, modify and distribute the contents of the project, given that proper references to the original author are provided.
