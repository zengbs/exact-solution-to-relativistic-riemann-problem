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


double FanFunction ( double Dens_at_Xi, void *params )
{
  struct Rarefaction *rarefaction = ( struct Rarefaction * ) params;
  
  double Xi        = rarefaction -> Xi;
  double Right_Yes = rarefaction -> Right_Yes;
  double sign = ( Right_Yes ) ? +1.0 : -1.0;

  double Velocity_at_Xi_1, Velocity_at_Xi_2;
  
  /* Step 1 */
  Velocity_at_Xi_1  = Isentropic_Dens2Velocity( Dens_at_Xi, params );

  /* Step 2 */
  double PresUp     = rarefaction -> PresUpStream;
  double DensUp     = rarefaction -> DensUpStream;
 
  double U_Xi       = Xi / sqrt( 1.0 - Xi*Xi );

  double Temp_at_Xi = Isentropic_Dens2Temperature( Dens_at_Xi, PresUp/DensUp, DensUp );

  double Cs         = Flu_SoundSpeed( Temp_at_Xi );

  RelativeVelocity( U_Xi, sign*Cs, NULL, &Velocity_at_Xi_2 );

  return Velocity_at_Xi_1 - Velocity_at_Xi_2;
}


double GetDensInFan( struct Rarefaction *Rarefaction )
{
 double Dens_at_Xi, DensUp, DensDown;

 DensUp   = Rarefaction -> DensUpStream;
 DensDown = Rarefaction -> DensDownStream;

 double DensMin = DensDown*0.9999;
 double DensMax = DensUp*1.0001;

 Dens_at_Xi = RootFinder( FanFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.5*(DensMax + DensMin), DensMin, DensMax );

 return Dens_at_Xi;
}

double GetPresInFan( double Dens_at_Xi, double PresUp, double DensUp )
{
  double Pres_at_Xi;

  Pres_at_Xi = Isentropic_Dens2Pres( Dens_at_Xi, PresUp/DensUp, DensUp );

  return Pres_at_Xi;
}

double GetVelocityInFan( double Xi,  double Dens_at_Xi, double Pres_at_Xi, bool Right_Yes )
{
  double gamma_Xi = 1.0 / sqrt( 1.0 - Xi*Xi );

  double U_Xi = Xi * gamma_Xi;
  
  double Velocity;

  double Tempertaure_at_Xi = Pres_at_Xi/Dens_at_Xi;

  double Cs = Flu_SoundSpeed( Tempertaure_at_Xi );

  if ( Right_Yes )
      RelativeVelocity( U_Xi, +Cs, NULL, &Velocity );
  else
      RelativeVelocity( U_Xi, -Cs, NULL, &Velocity );

  return Velocity;
}


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

  Temperature = RootFinder( Isentropic_TemperatureFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 5.0, 1.0, 10.0 );
  
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
         Cs      /= sqrt( 1.0 + Cs*Cs ); // 4-sound speed -> 3-sound speed

  double LorentzFactor = sqrt( 1.0 + y[0]*y[0] );
  f[0] = sign * LorentzFactor*Cs/Dens;

  return GSL_SUCCESS;
}



// density to the 4-velocity of flow

double Isentropic_Dens2Velocity ( double DensDown, struct Rarefaction *upstream )
{
  double VelyUp = upstream -> VelyUpStream;
  double DensUp = upstream -> DensUpStream;

  double t0 = DensUp;
  double t1 = DensDown;

  double ini_step = 1e-10;
  double abserr   = 0.0;
  double relerr   = __DBL_EPSILON__;

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
