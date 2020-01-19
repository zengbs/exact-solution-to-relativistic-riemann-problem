#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "Global.h"
#include "Prototypes.h"
#include "Macro.h"

void GetHeadTailVelocity( double PresUp, double DensUp, double VelocityUp,
			              double PresDown, double DensDown, double VelocityDown,
                          double *HeadVelocity, double *TailVelocity, bool Right_Yes )
{
  double Cs_Up, Cs_Down;

  Cs_Up   = Flu_SoundSpeed( PresUp / DensUp   ); 
  Cs_Down = Flu_SoundSpeed( PresDown / DensDown ); 


  if ( Right_Yes )
  {
	 *HeadVelocity = Cs_Up   * sqrt(1.0+SQR(VelocityUp))   + VelocityUp   * sqrt(1.0+SQR(Cs_Up));
	 *TailVelocity = Cs_Down * sqrt(1.0+SQR(VelocityDown)) + VelocityDown * sqrt(1.0+SQR(Cs_Down));
  }
  else
  {
	 *HeadVelocity = - Cs_Up   * sqrt(1.0+SQR(VelocityUp))   + VelocityUp   * sqrt(1.0+SQR(Cs_Up));
	 *TailVelocity = - Cs_Down * sqrt(1.0+SQR(VelocityDown)) + VelocityDown * sqrt(1.0+SQR(Cs_Down));
  }
}

double GetDensDownRarefaction( double PresDown, double PresUp, double DensUp )
{
  double DensDown;

  DensDown = pow( DensUp, Gamma )*PresDown / PresUp;

  DensDown = pow( DensDown, 1.0 / Gamma );
 
  return DensDown;
}


double GetVelocityDownRarefaction( double PresDown, double DensDown, double PresUp, double DensUp, double VelocityUp, bool Right_Yes )
{
  double Velocity;

  if ( Right_Yes )
  {
	double A_Minus;

	A_Minus  = A_MinusFun( PresUp / DensUp );
	A_Minus /= A_MinusFun( PresDown / DensDown);

    Velocity  = VelocityUp * ( A_Minus + 1.0 ) + sqrt( 1.0 + SQR(VelocityUp) ) * ( A_Minus - 1.0 );
	Velocity /= sqrt( 4.0 * A_Minus );
    Velocity *= SIGN( VelocityUp * ( A_Minus - 1.0 ) + sqrt( 1.0 + SQR(VelocityUp) ) * ( A_Minus + 1.0 ) );
  }
  else
  {
	double A_Plus;

	A_Plus  = A_PlusFun( PresUp / DensUp );
	A_Plus /= A_PlusFun( PresDown / DensDown);

    Velocity  = VelocityUp * ( A_Plus + 1.0 ) + sqrt( 1.0 + SQR(VelocityUp) ) * ( A_Plus - 1.0 );
	Velocity /= sqrt( 4.0 * A_Plus );
    Velocity *= SIGN( VelocityUp * ( A_Plus - 1.0 ) + sqrt( 1.0 + SQR(VelocityUp) ) * ( A_Plus + 1.0 ) );
  }

  return Velocity;
}

double GetDensInFan( double Cs, double PresUp, double DensUp )
{
  double k = PresUp * pow (DensUp, -Gamma);
  double tmp = ( 1.0 / (Cs * Cs)) - (1.0 / Gamma_1 );
  double Dens = pow ( k * Gamma * tmp, -1.0 / Gamma_1 );

  return Dens;
}

double GetPresInFan( double DensInFan, double PresUp, double DensUp )
{
  double Pres = PresUp * pow (DensInFan / DensUp, Gamma);

  return Pres;
}

double GetVelocityInFan( double Cs, double Xi, bool Right_Yes )
{
  double gamma_Xi = 1.0 / sqrt( 1.0 - Xi*Xi );

  double U_Xi = Xi * gamma_Xi;
  
  double Velocity;

  if ( Right_Yes )
	  Velocity = - Cs * gamma_Xi + sqrt(1.0+Cs*Cs) * U_Xi;
  else
	  Velocity = + Cs * gamma_Xi + sqrt(1.0+Cs*Cs) * U_Xi;

  return Velocity;
}


double GetSoundSpeedInFan ( struct Rarefaction *Rarefaction )
{
  double Temp, Cs;

  Temp = RootFinder( TemperatureFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 0.0, sqrt (Gamma_1) - 1e-16 );

  Cs = Flu_SoundSpeed( Temp );

  return Cs;
}

double TemperatureFunction ( double Temp, void *params )
{
  struct Rarefaction *Fan = ( struct Rarefaction * ) params;

  bool   Right_Yes    = Fan -> Right_Yes   ;
  double PresUp       = Fan -> PresUpStream;
  double DensUp       = Fan -> DensUpStream;
  double VelocityUp   = Fan -> VelyUpStream;
  double Xi           = Fan -> Xi          ;
 
  double Velocity, Var0, Var1, Cs_Up, Cs;

  double Sqrt_Gamma_1 = sqrt(Gamma_1);

  Cs_Up = Flu_SoundSpeed( PresUp / DensUp );

  Cs    = Flu_SoundSpeed( Temp );

  double gamma_Xi = 1.0 / sqrt( 1.0 - Xi*Xi );

  double U_Xi = Xi * gamma_Xi;

  if ( Right_Yes )
  {
    //Velocity = ( Xi - Cs )/( 1.0 - Cs * Xi );

    //Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	//Var0  = pow( Var0, -2.0 / Sqrt_Gamma_1 );

	//Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity );


	//Var1  = ( Sqrt_Gamma_1 + Cs_Up ) / ( Sqrt_Gamma_1 - Cs_Up );

	//Var1  = pow( Var1, -2.0 / Sqrt_Gamma_1 );

	//Var1 *= ( 1.0 + VelocityUp ) / ( 1.0 - VelocityUp );
	

    Velocity = - Cs * gamma_Xi + sqrt(1.0 + Cs*Cs) * U_Xi;

    double gamma_Velocity = sqrt( 1.0 + Velocity*Velocity );

	Var0  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );

	Var0  = pow( Var0, 4.0 / Sqrt_Gamma_1 );

	Var0 *= SQR( gamma_Velocity + Velocity );

    double gamma_VelocityUp = sqrt( 1.0 + VelocityUp*VelocityUp );

	Var1  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );

	Var1  = pow( Var1, 4.0 / Sqrt_Gamma_1 );

	Var1 *= SQR( gamma_VelocityUp + VelocityUp );
  }
  else
  {
    //Velocity = ( Xi + Cs )/( 1.0 + Cs * Xi );

    //Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	//Var0  = pow( Var0, +2.0 / Sqrt_Gamma_1 );

	//Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity );


	//Var1  = ( Sqrt_Gamma_1 + Cs_Up ) / ( Sqrt_Gamma_1 - Cs_Up );

	//Var1  = pow( Var1, +2.0 / Sqrt_Gamma_1 );

	//Var1 *= ( 1.0 + VelocityUp ) / ( 1.0 - VelocityUp );

    Velocity = + Cs * gamma_Xi + sqrt(1.0 + Cs*Cs) * U_Xi;

    double gamma_Velocity = sqrt( 1.0 + Velocity*Velocity );

	Var0  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );

	Var0  = pow( Var0, -4.0 / Sqrt_Gamma_1 );

	Var0 *= SQR( gamma_Velocity + Velocity );

    double gamma_VelocityUp = sqrt( 1.0 + VelocityUp*VelocityUp );

	Var1  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );

	Var1  = pow( Var1, -4.0 / Sqrt_Gamma_1 );

	Var1 *= SQR( gamma_VelocityUp + VelocityUp );
  }
  return Var1 - Var0;

}

