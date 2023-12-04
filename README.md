# Exact Solution to 1D Relativistic Riemann Problem

## Description

A tool for solving the exact solution to one-dimensional special relativistic Riemann problem.

## The avalibilities of Equation of State (EoS):

* Constant- $\Gamma$:

$$
\frac{h_{\Gamma}}{c^2}=1+\frac{\Gamma}{\Gamma-1}\left(\frac{k_{B}T}{mc^2}\right)
$$

* Taub-Mathews:

$$
\frac{h_{\text{TM}}}{c^2}=2.5\left(\frac{k_{B}T}{mc^2}\right)+\sqrt{2.25\left(\frac{k_{B}T}{mc^2}\right)^2+1}
$$

## Prerequisite Library
GNU Scientific Library (GSL)

## Quick Start
1. `git clone git@github.com:zengbs/exact-solution-to-relativistic-riemann-problem.git`
2. `cd exact-solution-to-relativistic-riemann-problem`
3. Edit the `GSL_DIR` in `Makefile` to match the path of GSL package
4. Edit `Input__Parameter` to set up initial conditions.
5. `make clean`
6. `make`
7. `./bin/a.out`
8. `cd plot_script`
9. `python plot__profiles.py`
10. `eog fig__profiles.png`

## Output Result
<img src="https://github.com/zengbs/exact-solution-to-relativistic-riemann-problem/blob/master/plot_script/fig__profiles.png" width="480">

## Reference
Please refer to Appendix C in [Tseng et al. 2021, MNRAS, 504, 3298](https://github.com/zengbs/published-papers/blob/main/2021-An_adaptive_mesh_GPU-accelerated_and_error_minimized_special_relativistic_hydrodynamics_code.pdf).
