#include <stdio.h>
#include "Prototypes.h"


int main()
{
  double DensLeft      =  1.0;
  double VelocityLeft  = -0.2;
  double PresLeft      =  2.0;

  double DensRight     =  1.0;
  double VelocityRight =  0.5;
  double PresRight     =  1.0

  struct InitialCondition IC = 
  {
     DensLeft     
     VelocityLeft 
     PresLeft     
                  
     DensRight    
     VelocityRight
     PresRight    
  };
  
  double PresStar; 
  
  PresStar = RootFinder( &IC, 0.0, __DBL_EPSILON__ );

  printf("PresStar=%e\n", PresStar);

  return 0;
}
