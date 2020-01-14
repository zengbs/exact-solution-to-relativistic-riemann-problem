#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"
#include "Macro.h"

void GetShockVelocity( double PresUp,   double DensUp,   double V_Up, 
		     		   double PresDown, double DensDown, double V_Down,
			           double *Vs_Left, double *Vs_Right )
{
  double ShockVelocity, LorentzFactor, J;
  double A, B, C;
  
  J = MassCurrent( PresUp, DensUp, PresDown, DensDown );
  
  if ( V_Up != V_Up   )
  {
    LorentzFactor = 1.0 / sqrt( 1.0 - SQR(V_Down) );

    A = SQR(J) + SQR( DensDown * LorentzFactor );
    B = -2.0 * V_Down * SQR( DensDown * LorentzFactor );
    C = SQR( DensDown*LorentzFactor*V_Down ) - SQR( J );
  }

  if ( V_Down != V_Down )
  {
    LorentzFactor = 1.0 / sqrt( 1.0 - SQR(V_Up) );

    A = SQR(J) + SQR( DensUp * LorentzFactor );
    B = -2.0 * V_Up * SQR( DensUp * LorentzFactor );
    C = SQR( DensUp*LorentzFactor*V_Up ) - SQR( J );
  }

  QuadraticSolver( A, B, C, Vs_Right, Vs_Left );
}

// solve eq. (4.140) for Va or Vb

double GetVelocityDown( double PresUp,   double DensUp, double ShockFrontVelocity,
                        double PresDown, double DensDown )
{
  double A, B, C, LorentzFactor, J, V_Left, V_Right;

  LorentzFactor = 1.0 / sqrt( 1.0 - SQR(ShockFrontVelocity) );

  J = MassCurrent( PresUp, DensUp, PresDown, DensDown );

  A = SQR( DensDown * LorentzFactor ) + SQR( J );

  B = -2.0 * ShockFrontVelocity * SQR( DensDown * LorentzFactor );

  C = SQR( DensDown * LorentzFactor * ShockFrontVelocity ) - SQR( J );

  QuadraticSolver( A, B, C, &V_Left, &V_Right );

  if( ShockFrontVelocity > 0.0 ) return V_Right; 
  else                           return  V_Left;

}


double GetDensDown( double PresUp, double DensUp, double PresDown  )
{
  double EnthalpyDown, DensDown;

  EnthalpyDown = TaubAdiabatic(PresUp, DensUp, PresDown);
  
  DensDown = ( Gamma / Gamma_1 ) * ( PresDown / ( EnthalpyDown - 1.0 ) );

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


double TaubAdiabatic ( double PresUp, double DensUp, double PresDown )
{
    double EnthalpyUp, EnthalpyDown, PresDiff;

	EnthalpyUp = Flu_Enthalpy( PresUp, DensUp );

    PresDiff = PresUp - PresDown;

    double A, B, C;

	A = 1.0 + Gamma_1 * PresUp / PresDown;
	B = - Gamma_1 * PresDiff / PresDown;
	C = - Gamma_1 * PresDiff * EnthalpyUp / PresUp - ( 1.0 + Gamma_1 * PresDown / PresUp )*SQR(EnthalpyUp);

    QuadraticSolver( A, B, C, &EnthalpyDown, NULL );

    return EnthalpyDown;
}
