# Exact Solution to 1D Relativistic Riemann Problem

Support Equation of State (EoS):
* Constant- $\Gamma$:

$$
\frac{h_{\Gamma}}{c^2}=1+\frac{\Gamma}{\Gamma-1}\left(\frac{k_{B}T}{mc^2}\right)
$$

* Taub-Mathews

$$
\frac{h_{\text{TM}}}{c^2}=2.5\left(\frac{k_{B}T}{mc^2}\right)+\sqrt{2.25\left(\frac{k_{B}T}{mc^2}\right)^2+1}
$$

## Prerequisite Library
GNU Scientific Library

## Quick Start
1. `git clone git@github.com:zengbs/exact-solution-to-relativistic-riemann-problem.git`
2. `cd exact-solution-to-relativistic-riemann-problem`
3. `make clean`
4. `make`
5. `./bin/a.out`
6. `cd plot`
7. `python plot__profiles.py`
8. `eog fig__profiles.png`

## Example
<img src="https://github.com/zengbs/exact-solution-to-relativistic-riemann-problem/blob/master/plot/fig__profiles.png" width="480">

## Reference
Please refer to Appendix C in the [code paper](https://github.com/zengbs/published-papers/blob/main/2021-An_adaptive_mesh_GPU-accelerated_and_error_minimized_special_relativistic_hydrodynamics_code.pdf).
