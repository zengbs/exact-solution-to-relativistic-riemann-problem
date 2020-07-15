#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"
#include "Macro.h"

void GetShockVelocity( double PresUp,   double DensUp,   double V_Up, 
		     		   double PresDown, double DensDown,
			           double *Vs_Left, double *Vs_Right )
{
  double J;
  
  J = MassCurrent( PresUp, DensUp, PresDown, DensDown );

  if ( Vs_Right != NULL )
     *Vs_Right = +( J * sqrt( 1.0 + V_Up*V_Up ) + V_Up * sqrt( J*J + DensUp*DensUp ) ) / DensUp;

  if ( Vs_Left != NULL )
     *Vs_Left  = -( J * sqrt( 1.0 + V_Up*V_Up ) - V_Up * sqrt( J*J + DensUp*DensUp ) ) / DensUp;
}


double GetVelocityDown( double PresUp,   double DensUp, double ShockVelocity,
                        double PresDown, double DensDown )
{
  double J;

  J = MassCurrent( PresUp, DensUp, PresDown, DensDown );

  if( ShockVelocity > 0.0 )
  {
     double Velocity_Left;

	 Velocity_Left = -J * sqrt( 1.0 + SQR(ShockVelocity) ) + ShockVelocity * sqrt( J*J + SQR(DensDown) );
	 Velocity_Left /= DensDown;

	 return Velocity_Left;
  }
  else
  {
     double Velocity_Right;

	 Velocity_Right = J * sqrt( 1.0 + SQR(ShockVelocity) ) + ShockVelocity * sqrt( J*J + SQR(DensDown) );
	 Velocity_Right /= DensDown;

	 return Velocity_Right;
  }

}


double GetDensDown( double PresUp, double DensUp, double PresDown  )
{
  double EnthalpyDown, DensDown;

  EnthalpyDown = GetEnthalpyDown(PresUp, DensUp, PresDown);
  double TempDown;
  TempDown = Enthalpy2Temperature( EnthalpyDown );
  DensDown = PresDown/TempDown;

  return DensDown;
}

double MassCurrent( double PresUp, double DensUp, double PresDown, double DensDown )
{
  double MassCurrent;
  double EnthalpyUp, EnthalpyDown;

  EnthalpyUp   = Flu_Enthalpy(   PresUp,   DensUp );
  EnthalpyDown = Flu_Enthalpy( PresDown, DensDown );
  MassCurrent  = PresDown - PresUp;
  MassCurrent /= ( EnthalpyUp / DensUp ) - ( EnthalpyDown / DensDown );

  if ( MassCurrent < 0.0 )
  {
    printf("MassCurrent is %e!!\n", MassCurrent);
    exit(1);
  }

  MassCurrent = sqrt( MassCurrent );
  
  return MassCurrent;
}

struct Parameters
{
  double EnthalpyUp;
  double PresUp    ;
  double DensUp    ;
  double PresDown  ;
};

double GetEnthalpyDown ( double PresUp, double DensUp, double PresDown )
{
    double EnthalpyUp, EnthalpyDown;

	EnthalpyUp = Flu_Enthalpy( PresUp, DensUp );
#   if ( EOS == GAMMA )
    double PresDiff = PresUp - PresDown;

    double A, B, C;

	A = 1.0 + Gamma_1 * PresUp / PresDown;
	B = - Gamma_1 * PresDiff / PresDown;
	C = - Gamma_1 * PresDiff * EnthalpyUp / PresUp - ( 1.0 + Gamma_1 * PresDown / PresUp )*SQR(EnthalpyUp);
 

    QuadraticSolver( A, B, C, &EnthalpyDown, NULL );
#   else

    struct Parameters params;

    params.EnthalpyUp = EnthalpyUp;
    params.PresUp     = PresUp;
    params.DensUp     = DensUp;
    params.PresDown   = PresDown;

    EnthalpyDown = RootFinder( JumpConditionForEnthalpy, (void*)&params, 0.0, __DBL_EPSILON__, 1e6, 1e5, 1e7, __FUNCTION__  );
#   endif
    return EnthalpyDown;
}

double JumpConditionForEnthalpy( double EnthalpyDown, void* params )
{
    struct Parameters *pparams = (struct Parameters *) params;
  
    double EnthalpyUp = pparams -> EnthalpyUp; 
    double PresUp     = pparams -> PresUp    ;
    double DensUp     = pparams -> DensUp    ;
    double PresDown   = pparams -> PresDown  ;
    double TempDown   = Enthalpy2Temperature( EnthalpyDown );
    double TempUp     = PresUp/DensUp;

    double LeftSide   = EnthalpyUp*EnthalpyUp - EnthalpyDown*EnthalpyDown;
    double RightSide  = ( TempUp - PresDown/DensUp )*EnthalpyUp + ( PresUp/PresDown - 1.0 )*EnthalpyDown*TempDown;
    return LeftSide - RightSide;
}
