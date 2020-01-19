#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"

//double Gamma = 5.0/3.0;
//double Gamma_1 = 0.666666666666666666666;
double Gamma = 4.0/3.0;
double Gamma_1 = 0.333333333333333333333;


int main()
{
  double DensLeft      =  0.2;
  double VelocityLeft  =  0.0;
  double PresLeft      =  0.5;

  double DensRight     =  1.0;
  double VelocityRight =  -1.0;
  double PresRight     =  1.0;

  
  double DT            = 0.1;
  double End_T         = 0.5;
  double X_Left        = 0.0;
  double X_Right       = 1.0;
  int NCell            = 1024;
 

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

  //double up, lb;

  //up = 1e2;
  //lb = 1e-2;
  //int N = 1000;

  //double dp = (up-lb)/(double)N;

  //double pres = 0.0;

  //double fun_pres;

  //for (int i=1;i<=N;i++)
  //{
  //  pres = i*dp;
  //  fun_pres = PresFunction( pres, &IC );
  //
  //  printf("%e  %e\n", pres, fun_pres);
  //}


  //double up, lb;
  //up = 1e13;
  //lb = 1e12;
  //int N = 1000;

  //double dp = (up-lb)/(double)N;

  //double pres = 0.0;

  //double fun_pres;

  //for (int i=1;i<=N;i++)
  //{
  //  pres = i*dp;
  //  fun_pres = PresFunction( pres, &IC );
  //
  //     printf("%e  %e\n", pres, fun_pres);
  //}



  return 0;
}
