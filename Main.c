#include <stdio.h>
#include <stdbool.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"

double Gamma = 5.0/3.0;
double Gamma_1 = 0.6666666666666666;

int main()
{

  double DensLeft      =  2.0;
  double VelocityLeft  =  0.0;
  double PresLeft      =  0.04;

  double DensRight     =  1.0;
  double VelocityRight =  0.0;
  double PresRight     =  0.02;

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
