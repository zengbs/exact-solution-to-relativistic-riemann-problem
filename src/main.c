#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include "../include/global.h"
#include "../include/prototypes.h"
#include "../include/struct.h"

// the value is given in Makefile
double Gamma = GAmma;
double Gamma_1 = GAmma_1;

// input parameters
double L_X, L_DENS, L_VELX, L_PRES;
double R_X, R_DENS, R_VELX, R_PRES;
double DT, END_T;
int N_CELL;

int main() {
  Load_Parameter();

  printf("===========================================================\n");
  printf("The input LR states are:\n");
  printf("L_X                 = %24.16e\n", L_X);
  printf("L_DENS              = %24.16e\n", L_DENS);
  printf("L_VELX              = %24.16e\n", L_VELX);
  printf("L_PRES              = %24.16e\n", L_PRES);
  printf("R_X                 = %24.16e\n", R_X);
  printf("R_DENS              = %24.16e\n", R_DENS);
  printf("R_VELX              = %24.16e\n", R_VELX);
  printf("R_PRES              = %24.16e\n", R_PRES);
  printf("===========================================================\n");
  printf("The other parameters are:\n");
  printf("DT                  = %24.16e\n", DT);
  printf("END_T               = %24.16e\n", END_T);
  printf("N_CELL              = %d\n", N_CELL);
  printf("===========================================================\n");

  struct InitialCondition IC = {
      L_DENS, L_VELX, L_PRES, R_DENS, R_VELX, R_PRES,
  };

  struct RiemannProblem RP;

  struct PlotParams plot = {
      DT, END_T, L_X, R_X, N_CELL,
  };

  int Pattern;

  Pattern = GetAllInfomation(&IC, &RP);

  Plot(Pattern, &RP, plot);

  return 0;
}  // FUNCTION : main
