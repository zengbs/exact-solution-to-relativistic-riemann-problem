# ExactSolutionRelativisticRiemannProblem
This is a tool of solving the exact solution of special relativistic Riemann problem.

Author: [Po-Hsun Tseng](https://github.com/zengbs)


## Usage
1. If necessary, please change the `GSL_DIR` in `Makefile` to match your machine setup.
   ```makefile
   # paths
   SRC_DIR := ./src
   OBJ_DIR := ./obj
   INC_DIR := ./include
   GSL_DIR := /software/gsl/default
   ```
1. Set the fluid initial condition at the begining of `int main()` in `./src/Main.c`. Please check out [Tseng et al. 2021, MNRAS, 504, 3298](https://academic.oup.com/mnras/article-abstract/504/3/3298/6224873?redirectedFrom=PDF) for more details.

   ```c
   double DensLeft       = 1.e-5;
   double VelocityLeft   = 1.e+6;
   double PresLeft       = 1.0;
   
   double DensRight      = 1.e-5;
   double VelocityRight  = -1.e+6;
   double PresRight      = 1.0;
   ```
1. Compile the software.
   ```shell
   > make clean
   > make
   ...
   ...
   Linking to SR_Riemann_solver ... Successful!
   ```
1. Execute it!
   ```shell
   > ./SR_Riemann_solver
   SS pattern !!
   PresStar=5.3333333333797783e+12
   Pattern=1
   ```
1. There should be `*.dat` in current directory after executing. You may use the simple script under `scripts/plot_dat.py` to plot the solution.
