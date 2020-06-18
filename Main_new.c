#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"


double Gamma = 5.0/3.0;
double Gamma_1 = 0.666666666666666666666;
//double Gamma = 4.0/3.0;
//double Gamma_1 = 0.333333333333333333333;


int main()
{
  double DensLeft      =  1.0;
  double VelocityLeft  =  0.0;
  double PresLeft      =  10.0;

  double DensRight     =  1.0;
  double VelocityRight = -1.0;
  double PresRight     =  1.0;
//// RR
//  double DensLeft      =  1.0;
//  double VelocityLeft  = -1.0;
//  double PresLeft      =  0.1;
//
//  double DensRight     =  1.0;
//  double VelocityRight = +2.0;
//  double PresRight     =  1.0;
//
//// RS
//  double DensLeft      =  1.0;
//  double VelocityLeft  =  0.0;
//  double PresLeft      =  1.0;
//
//  double DensRight     =  1.0;
//  double VelocityRight =  0.0;
//  double PresRight     =  0.1;
//
//// SR
//  double DensLeft      =  1.0;
//  double VelocityLeft  =  0.0;
//  double PresLeft      =  0.1;
//
//  double DensRight     =  1.0;
//  double VelocityRight =  0.0;
//  double PresRight     =  1.0;
//
////SS
//  double DensLeft      =  5.0;
//  double VelocityLeft  = +1.0;
//  double PresLeft      = 10.0;
//
//  double DensRight     = 10.0;
//  double VelocityRight = -1.0;
//  double PresRight     =  1.0;



  
  double DT            = 0.1;
  double End_T         = 0.6;
  double X_Left        = 0.0;
  double X_Right       = 1.0;
  int NCell            = 4096;
 

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
