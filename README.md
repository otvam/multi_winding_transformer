#  Optimal Current Distribution in Multi-Winding Transformers

## Introduction

This repository contains an implementation of the **optimal current sharing** in **multi-port transformers**.
From the **impedance matrix**, the **optimal currents** (amplitudes and phases) is computed in order to **minimize the losses**.

A **full description** of the optimal operation of **multi-port transformers** can be found in:

* **Optimal Current Distribution in Multi-Winding Transformers for Isolated and Wireless Power Transfer**
* **K. Datta, Y. Wu, C. R. Sullivan, and J. T. Stauth**
* **https://doi.org/10.1109/OJPEL.2025.3590020**
* **IEEE OJPEL 2025**

Two different **optimization methods** are available:
* **Eigenvalue method** (neglect the mutual resistances, described in the paper)
* **Numerical solver** (include the mutual resistances, not part of the paper)

The **numerical solver** features several advantages compared to the **eigenvalue method**:
* Usage of a standard quadratically constrained quadratic solver.
* Analytical expression for the gradient vectors and hessian matrices.
* Complex objective function including the reactive power.
* Inclusion of the mutual resistance coefficients.

## Main Files

* [run_single.m](run_single.m) - Solve the optimal current sharing problem for a multi-winding transformer 
* [run_sweep.m](run_sweep.m) - Optimal current sharing with loss and reactive power minimization
* [get_problem.m](get_problem.m) - Definition of the transformer and the power flow
* [get_tolerance.m](get_tolerance.m) - Definition of the numerical tolerance

## Compatibility

* Tested with MATLAB R2024b.
* The `gads_toolbox` is required.
* The `optimization_toolbox` is required.

## Author

* **Thomas Guillod**
* Email: guillod@otvam.ch
* Website: https://otvam.ch

## Credits

This code was created at **Dartmouth College** by the research group of **Prof. Sullivan**:
* Dartmouth College, NH, USA: https://dartmouth.edu
* Dartmouth Engineering: https://engineering.dartmouth.edu
* PMIC: https://pmic.engineering.dartmouth.edu

## License

This project is licensed under the **MIT License**, see [LICENSE.md](LICENSE.md).
