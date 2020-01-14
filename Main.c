#include <stdio.h>
#include <stdbool.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"

double Gamma = 5.0/3.0;
double Gamma_1 = 0.666666666666666666666;
//double Gamma = 4.0/3.0;
//double Gamma_1 = 0.333333333333333333333;


int main()
{
  double DensLeft      =  0.5;
  double VelocityLeft  =  4.47213595499957928e-01;
  double PresLeft      =  0.04;

  double DensRight     =  0.5;
  double VelocityRight = -4.47213595499957928e-01;
  double PresRight     =  0.04;

  
  double DT            = 0.1;
  double End_T         = 0.5;
  double X_Left        = 0.0;
  double X_Right       = 1.0;
  int NCell            = 512;
 


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

  Plot( Pattern, &RP, plot );

  return 0;
}
