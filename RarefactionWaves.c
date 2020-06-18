#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "Global.h"
#include "Prototypes.h"
#include "Macro.h"


static double Isentropic_TemperatureFunction ( double Temperature, void *params );

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

  Temp = RootFinder( TemperatureFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 1e-3, 6.0 );

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
 
  double Velocity, Var0, Var1, Cs, TempUp;

  double Sqrt_Gamma_1 = sqrt(Gamma_1);

  TempUp = PresUp / DensUp;

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
// ============================================= OK
double Isentropic_Constant ( double Init_Temp, double Init_Dens )
{
  double K;
# if ( EOS == GAMMA )
  K = Init_Temp / pow( Init_Dens, Gamma_1 );
# elif ( EOS == TM )
  K  = Init_Temp*( 1.5*Init_Temp + sqrt( 2.25*Init_Temp*Init_Temp + 1.0 ) );
  
  K /= pow(Init_Dens, 2.0/3.0);
# endif
  return K;
}

//========================================= ??

double Isentropic_Dens2Temperature ( double Dens, double Init_Temp, double Init_Dens )
{
  double Temperature, K;
  K = Isentropic_Constant(Init_Temp, Init_Dens);

# if ( EOS == GAMMA )
 Temperature = K*pow( Dens, Gamma_1 );
# elif ( EOS == TM )
  double A = K*pow(Init_Dens, 2.0/3.0);
  Temperature = A / sqrt( 3.0*A + 1.0 );
# endif

  return Temperature;
}

double Isentropic_Temperature2Dens ( double Temperature, double Init_Temp, double Init_Dens )
{
  double Dens, K;
  K = Isentropic_Constant(Init_Temp, Init_Dens);

# if ( EOS == GAMMA )
  Dens = pow( Temperature/K, 1.0/Gamma_1 );

# elif ( EOS == TM )
  Dens  = 1.5*Temperature*Temperature + Temperature*sqrt(2.25*Temperature*Temperature + 1.0);

  Dens /= K;

  Dens  = pow( Dens, 1.5 );
# endif

  return Dens;
}

//========================================= OK
double Isentropic_Pres2Temperature ( struct Rarefaction *Rarefaction )
{
  double Temperature;

  Temperature = RootFinder( Isentropic_TemperatureFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 1e-2, 10.0 );
  
  return Temperature;
}


//
// Pres = Pres( Temp ) OK
//
double Isentropic_TemperatureFunction ( double TempDown, void *params ) 
{
  struct Rarefaction *rarefaction = ( struct Rarefaction * ) params;

  double PresUp    = rarefaction -> PresUpStream;
  double DensUp    = rarefaction -> DensUpStream;
  double PresDown  = rarefaction -> PresDownStream;
  double TempUp    = PresUp / DensUp;
  double K         = Isentropic_Constant(TempUp, DensUp);
  double Expression;
//  printf("K=%e\n", K);
# if ( EOS == GAMMA )
  Expression  = pow( TempDown, Gamma/Gamma_1 );
  Expression *= pow( K, -1.0/Gamma_1 );
  
# elif ( EOS == TM )
  Expression  = pow( TempDown, 5.0/3.0 )*( 1.5*TempDown + sqrt(2.25*TempDown*TempDown + 1.0) );

  Expression /= K;
  Expression *= pow( Expression, 1.5 );
# endif

  return Expression - PresDown;
}

double Isentropic_Temperature2Pres ( double Temperature, void *params  )
{
  struct Rarefaction *upstream = ( struct Rarefaction * ) params;

  double Init_Pres  = upstream -> PresUpStream;
  double Init_Dens  = upstream -> DensUpStream;
  double Init_Temp  = Init_Pres / Init_Dens;
  double Pres, K;

  K = Isentropic_Constant(Init_Temp, Init_Dens);

# if ( EOS == GAMMA )
  Pres = pow( Temperature, Gamma )/K;
  Pres = pow( Pres, 1.0/Gamma_1 );
  
# elif ( EOS == TM )
  Pres  = pow( Temperature, 5.0/3.0 )*( 1.5*Temperature + sqrt(2.25*Temperature*Temperature + 1.0) );

  Pres /= K;

  Pres  = pow( Pres, 1.5 );
# endif

  return Pres;
}

//========================================= OK

double Isentropic_Pres2Dens ( struct Rarefaction *Rarefaction )
{
  double Temperature = Isentropic_Pres2Temperature( Rarefaction );
  double Pres        = Rarefaction -> PresDownStream;
  return Pres / Temperature;
}


double Isentropic_Dens2Pres ( double Dens, double Init_Temp, double Init_Dens )
{
  double Pres;

  Pres = Dens * Isentropic_Dens2Temperature( Dens, Init_Temp, Init_Dens );

  return Pres;
}


//============ Solve ODE =============================

// dU/d rho = LorentzFactor * Cs / rho
int func ( double Dens, const double y[], double f[], void *params )
{
  struct Rarefaction *upstream = (struct Rarefaction *)params;

  bool Right_Yes = upstream->Right_Yes;
  double sign = ( Right_Yes ) ? +1.0 : -1.0;

  double DensUp   = upstream->DensUpStream;
  double PresUp   = upstream->PresUpStream;
  double TempUp   = PresUp / DensUp;

  double TempDown = Isentropic_Dens2Temperature( Dens, TempUp, DensUp );
  double Cs       = Flu_SoundSpeed( TempDown );

  double LorentzFactor = sqrt( 1.0 + y[0]*y[0] );
  f[0] = sign * LorentzFactor*Cs/Dens;

  return GSL_SUCCESS;
}



double Isentropic_Dens2Velocity ( double DensDown, struct Rarefaction *upstream )
{
  double VelyUp = upstream -> VelyUpStream;
  double DensUp = upstream -> DensUpStream;

  double t0 = DensUp;
  double t1 = DensDown;

  double ini_step = 1e-10;
  double abserr   = 0.0;
  double relerr   = 1e-16;

  if ( t1 < t0 ) ini_step *= -1.0;

  gsl_odeiv2_system sys = {func, NULL, 1, upstream};

  gsl_odeiv2_driver * d =  gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, ini_step, abserr, relerr);

  double y[1] = { VelyUp };

  int status = gsl_odeiv2_driver_apply (d, &t0, t1, y);

  if (status != GSL_SUCCESS)
  {
      printf ("error, return value=%d\n", status);
      exit(0);
  }

  gsl_odeiv2_driver_free (d);

  return y[0]; // return VelyDown;
}
