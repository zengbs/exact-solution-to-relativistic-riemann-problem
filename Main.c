#include <stdio.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"

double Gamma = 1.33333333333333;
double Gamma_1 = Gamma - 1.0;

int main()
{

  double DensLeft      =  1.0;
  double VelocityLeft  = -0.2;
  double PresLeft      =  2.0;

  double DensRight     =  1.0;
  double VelocityRight =  0.5;
  double PresRight     =  1.0;

  struct InitialCondition IC = 
  {
     DensLeft,
     VelocityLeft,
     PresLeft,
     DensRight,    
     VelocityRight,
     PresRight,
  };
  
  double PresStar; 
  
  PresStar = RootFinder( &IC, 0.0, __DBL_EPSILON__ );

  printf("PresStar=%e\n", PresStar);

  return 0;
}
