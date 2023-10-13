#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"


double Gamma = GAmma;
double Gamma_1 = GAmma_1;


int main()
{
   double DensLeft       = 1.e-5;
   double VelocityLeft   = 1.e+6;
   double PresLeft       = 1.0;

   double DensRight      = 1.e-5;
   double VelocityRight  = -1.e+6;
   double PresRight      = 1.0;

   double DT            = 1.0;
   double End_T         = 1.0;
   double X_Left        = 0.0;
   double X_Right       = 1.0;
   int    NCell         = 4096;

   struct InitialCondition IC =
   {
      DensLeft,
      VelocityLeft,
      PresLeft,
      DensRight,
      VelocityRight,
      PresRight,
   };


   struct RiemannProblem RP;

   struct PlotParams plot =
   {
      DT,
      End_T,
      X_Left, X_Right,
      NCell,
   };

   int Pattern;

   Pattern = GetAllInfomation( &IC, &RP );

   printf("Pattern=%d\n", Pattern);
   Plot( Pattern, &RP, plot );
   return 0;
}
