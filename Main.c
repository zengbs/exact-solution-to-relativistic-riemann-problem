#include <stdio.h>
#include <stdbool.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"

double Gamma = 5.0/3.0;
double Gamma_1 = 0.6666666666666666;

int main()
{

  double DensLeft      =  0.5;
  double VelocityLeft  =  0.6;
  double PresLeft      =  0.02;

  double DensRight     =  2.0;
  double VelocityRight = -0.5;
  double PresRight     =  0.04;

 printf("%e\n", Flu_Enthalpy( PresLeft,  DensLeft  ) - PresLeft/DensLeft);
 printf("%e\n", Flu_Enthalpy( PresRight, DensRight ) - PresRight/DensRight);

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

  int Pattern;


  Pattern = GetAllInfomation( &IC, &RP );

  printf("Pattern=%d\n", Pattern);

  printf("%e\n", RP.SS.Leftt.ShockVelocity  );
  printf("%e\n", RP.SS.Leftt.PresUpStream   ); 
  printf("%e\n", RP.SS.Leftt.DensUpStream   ); 
  printf("%e\n", RP.SS.Leftt.VelyUpStream   ); 
  printf("%e\n", RP.SS.Leftt.PresDownStream ); 
  printf("%e\n", RP.SS.Leftt.DensDownStream ); 
  printf("%e\n", RP.SS.Leftt.VelyDownStream ); 


  return 0;
}
