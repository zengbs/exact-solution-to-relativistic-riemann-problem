# ExactSolutionRelativisticRiemannProblem
This is a tool for solving the exact solution of a special relativistic Riemann problem.

Author: [Po-Hsun Tseng](https://github.com/zengbs)


## Usage
1. Please change the `GSL_DIR` in `Makefile` to match your machine setup if necessary.
   ```makefile
   # paths
   SRC_DIR := ./src
   OBJ_DIR := ./obj
   INC_DIR := ./include
   GSL_DIR := /software/gsl/default
   ```
1. Compile the software.
   ```shell
   > make clean
   > make
   ...
   ...
   Linking to SR_Riemann_solver ... Successful!
   ```
1. Edit `Input__Parameter` to set the fluid initial condition. Please check out [Tseng et al. 2021, MNRAS, 504, 3298](https://academic.oup.com/mnras/article-abstract/504/3/3298/6224873?redirectedFrom=PDF) for more details.

   ```shell
   L_X       0.0      # left edge
   L_DENS    1.e-5    # density on the left state
   L_VELX    1.e+6    # velocity on the left state
   L_PRES    1.0      # pressure on the left state

   R_X       1.0      # right edge
   R_DENS    1.e-5    # density on the right state
   R_VELX    -1.e+6   # velocity on the right state
   R_PRES    1.0      # pressure on the right state

   DT        1.0      # physical output time increment
   END_T     1.0      # physical end time
   N_CELL    4096
   ```
1. Execute it!
   ```shell
   > ./SR_Riemann_solver
   ===========================================================
   The input LR states are:
   L_X                 =   0.0000000000000000e+00
   L_DENS              =   1.0000000000000001e-05
   L_VELX              =   1.0000000000000000e+06
   L_PRES              =   1.0000000000000000e+00
   R_X                 =   1.0000000000000000e+00
   R_DENS              =   1.0000000000000001e-05
   R_VELX              =  -1.0000000000000000e+06
   R_PRES              =   1.0000000000000000e+00
   ===========================================================
   The other parameters are:
   DT                  =   1.0000000000000000e+00
   END_T               =   1.0000000000000000e+00
   N_CELL              = 4096
   ===========================================================
   Shock-Shock pattern !!
   PresStar=5.3333333333797783e+12
   Plotting data 000000 (t=  0.0000000000000000e+00) ... Done
   ...
   ```
1. The current directory should have `*.dat` after executing. To plot the solution, you may use the simple script under `scripts/plot_dat.py`.
