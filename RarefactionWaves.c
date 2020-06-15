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
  double tmp;

  tmp  = Gamma_1 - Cs*Cs*(2.0-Gamma);
  tmp /= Gamma_1 + Cs*Cs*Gamma_1;
  tmp *= 1.0 + Cs*Cs;
  tmp /= Cs*Cs;

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

  Temp = RootFinder( TemperatureFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 1e-5, 1e2 );

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
 
  double Velocity, Var0, Var1, Cs_Up, Cs, TempUp;

  double Sqrt_Gamma_1 = sqrt(Gamma_1);

  TempUp = PresUp / DensUp;

  Cs_Up = Flu_SoundSpeed( TempUp );

  Cs    = Flu_SoundSpeed( Temp );

  double gamma_Xi = 1.0 / sqrt( 1.0 - Xi*Xi );

  double U_Xi = Xi * gamma_Xi;

  if ( Right_Yes )
  {
    Velocity = - Cs * gamma_Xi + sqrt(1.0 + Cs*Cs) * U_Xi;

    double gamma_Velocity = sqrt( 1.0 + Velocity*Velocity );

	Var0  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );

	Var0  = pow( Var0, 4.0 / Sqrt_Gamma_1 );

	Var0 *= SQR( gamma_Velocity + Velocity );

    double gamma_VelocityUp = sqrt( 1.0 + VelocityUp*VelocityUp );

	Var1  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * TempUp ) + sqrt( Gamma * TempUp ) );

	Var1  = pow( Var1, 4.0 / Sqrt_Gamma_1 );

	Var1 *= SQR( gamma_VelocityUp + VelocityUp );
  }
  else
  {
    Velocity = + Cs * gamma_Xi + sqrt(1.0 + Cs*Cs) * U_Xi;

    double gamma_Velocity = sqrt( 1.0 + Velocity*Velocity );

	Var0  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );

	Var0  = pow( Var0, -4.0 / Sqrt_Gamma_1 );

	Var0 *= SQR( gamma_Velocity + Velocity );

    double gamma_VelocityUp = sqrt( 1.0 + VelocityUp*VelocityUp );

	Var1  = Sqrt_Gamma_1 / ( sqrt( Gamma_1 + Gamma * TempUp ) + sqrt( Gamma * TempUp ) );

	Var1  = pow( Var1, -4.0 / Sqrt_Gamma_1 );

	Var1 *= SQR( gamma_VelocityUp + VelocityUp );
  }
  return Var1 - Var0;

}

double Isentropic_Constant ( double Init_Temp, double Init_Dens )
{
  double K;

  K  = Init_Temp*( 1.5*Init_Temp + sqrt( 2.25*Init_Temp*Init_Temp + 1.0 ) );
  
  K /= pow(Init_Dens, 2.0/3.0);

  return K;
}

double Isentropic_Dens2Temperature ( double Dens, double Init_Temp, double Init_Dens )
{
  double Temperature, A;

# if ( EOS == GAMMA )
# elif ( EOS == TM )
  A = pow(Init_Dens, 2.0/3.0) * Isentropic_Constant(Init_Temp, Init_Dens);
  Temperature = A / sqrt( 3.0*A + 1.0 );
# endif  
  
  return Temperature;
}

double Isentropic_Temperature2Pres ( double Dens, double Init_Temp, double Init_Dens )
{
  double Pres;

  Pres = Dens * Isentropic_Dens2Temperature( Dens, Init_Temp, Init_Dens );

  return Pres;
}

double Isentropic_Pres2Temperature ( double Pres, double Init_Temp, double Init_Dens )
{
  Temperature = RootFinder( Isentropic_Temperature2Pres, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 1e-5, 1e2 );

  return Temperature;
}

double Isentropic_Pres2Dens ( double Pres )
{
  double Temperature = Isentropic_Pres2Temperature( Pres, Init_Temp, Init_Dens );

  return Pres / Temperature;
}

// dU/d rho = LorentzFactor * Cs / rho
int RiemannInvariant ( double Dens, const double y[], double f[], void *params )
{
  double LorentzFactor = sqrt( 1.0 + Velocity*Velocity );

  double Temperature = Isentropic_Dens2Temperature( Dens, Init_Temp, Init_Dens );

  double Cs = Flu_SoundSpeed( Temperature );

  f[0] = LorentzFactor*Cs/Dens;

  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

