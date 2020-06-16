#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
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

//=========================================

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

//=========================================
double Isentropic_Pres2Temperature ( double Pres, struct Rarefaction *Rarefaction )
{
  double Temperature;

  Temperature = RootFinder( Isentropic_Temperature2Pres, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 1e-5, 1e2 );

  return Temperature;
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

//=========================================

double Isentropic_Pres2Dens ( double Pres, struct Rarefaction *Rarefaction )
{
  double Temperature = Isentropic_Pres2Temperature( Pres, Rarefaction );

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

  double LorentzFactor = sqrt( 1.0 + y[0]*y[0] );

  double Init_Dens = upstream->DensUpStream;
  double Init_Pres = upstream->PresUpStream;
  double Init_Temp = Init_Pres / Init_Dens;

  double Temperature = Isentropic_Dens2Temperature( Dens, Init_Temp, Init_Dens );

  double Cs = Flu_SoundSpeed( Temperature );

  f[0] = sign * LorentzFactor*Cs/Dens;

  return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  struct Rarefaction *upstream = (struct Rarefaction *)params;

  bool Right_Yes = upstream->    Right_Yes;
  double  DensUp = upstream-> DensUpStream;
  double  PresUp = upstream-> PresUpStream;
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 1, 1);

  gsl_matrix * m = &dfdy_mat.matrix;

  double sign;

  sign = ( Right_Yes ) ? +1.0 : -1.0;

  double LorentzFactor = sqrt(1.0 + y[0]*y[0]);

  double Temperature = Isentropic_Dens2Temperature( t, PresUp/DensUp, DensUp );

  double Cs = Flu_SoundSpeed( Temperature );

//  ∂ ⎛Cₛ⋅γ⎞        Cₛ⋅U    
//  ──⎜────⎟ = ─────────────
//  ∂U⎝ ρ  ⎠        ________
//                 ╱  2     
//             ρ⋅╲╱  U  + 1 
  double Jacobian = sign * ( y[0]/LorentzFactor ) * ( Cs / t );

  gsl_matrix_set (m, 0, 0, Jacobian);

  double K;
  K    = Isentropic_Constant( PresUp/DensUp, DensUp );


# if ( EOS == GAMMA )
//                     __________________                                
//                    ╱       Γ                                          
//                   ╱   Γ⋅K⋅ρ ⋅(Γ - 1)   ⎛ 2            Γ              ⎞
//             γ⋅   ╱   ──────────────── ⋅⎝Γ ⋅ρ - 2⋅Γ⋅K⋅ρ  - 4⋅Γ⋅ρ + 3⋅ρ⎠
//                 ╱         Γ                                           
//  ∂ ⎛Cₛ⋅γ⎞     ╲╱     Γ⋅K⋅ρ  + Γ⋅ρ - ρ                                 
//  ──⎜────⎟ = ──────────────────────────────────────────────────────────
//  ∂ρ⎝ ρ  ⎠                       2 ⎛     Γ          ⎞                  
//                              2⋅ρ ⋅⎝Γ⋅K⋅ρ  + Γ⋅ρ - ρ⎠                  

  dfdt[0]  = sqrt( Gamma * (Gamma - 1.0) * K * pow( t, Gamma ) );
  dfdt[0] /= sqrt( Gamma * K * pow(t, Gamma) + Gamma*t - t );
  dfdt[0] *= Gamma*Gamma*t - 2.0*Gamma*K*pow(t, Gamma) - 4.0*Gamma*t + 3.0*t;
  dfdt[0] /= 2.0*t*t*( Gamma*K*pow(t, Gamma) + Gamma*t - t );
# elif ( EOS == TM )
  double  dCs_dT, dT_dA, dA_drho, sqrtTemp, A, Temp;
  Temp = Isentropic_Dens2Temperature( t, PresUp/DensUp, DensUp );
  sqrtTemp = sqrt( 2.25*Temp*Temp + 1.0 );
  A        = pow(t, 2.0/3.0) * K;

//                    _______________________________________                                               
//                   ╱       ⎛             _____________⎞                                                   
//                  ╱        ⎜            ╱       2     ⎟     ⎛                       _____________        ⎞
//                 ╱       T⋅⎝4.5⋅T + 5⋅╲╱  2.25⋅T  + 1 ⎠     ⎜        2             ╱       2             ⎟
//                ╱    ───────────────────────────────────── ⋅⎝70.875⋅T  + 60.75⋅T⋅╲╱  2.25⋅T  + 1  + 33.75⎠
//               ╱                         _____________                                                    
//              ╱            2            ╱       2                                                         
//   d        ╲╱       18.0⋅T  + 12.0⋅T⋅╲╱  2.25⋅T  + 1  + 3                                                
//   ──(Cₛ) = ──────────────────────────────────────────────────────────────────────────────────────────────
//   dT           ⎛                        _____________                            _____________       ⎞   
//                ⎜        4          3   ╱       2                 2              ╱       2            ⎟   
//              T⋅⎝1458.0⋅T  + 972.0⋅T ⋅╲╱  2.25⋅T  + 1  + 799.875⋅T  + 330.75⋅T⋅╲╱  2.25⋅T  + 1  + 67.5⎠   

  dCs_dT   = 70.875*Temp*Temp + 60.75*Temp*sqrtTemp + 33.75;
  dCs_dT  *= sqrt(  4.5*Temp*Temp +  5.0*Temp*sqrtTemp );
  dCs_dT  /= sqrt( 18.0*Temp*Temp + 12.0*Temp*sqrtTemp + 3.0);
  dCs_dT  /= 1458.0*Temp*Temp*Temp*Temp + 972.0*Temp*Temp*Temp*sqrtTemp + 799.875*Temp*Temp + 330.75*Temp*sqrtTemp + 67.5;
  dCs_dT  /= Temp;

//   d          3⋅A + 2    
//   ──(T) = ──────────────
//   dA                 3/2
//           2⋅(3⋅A + 1)   

  dT_dA    = 3.0*A+2.0;
  dT_dA   /= 2.0*pow(3.0*A + 1.0 , 1.5);

//   d                            -0.333333333333333
//   ──(A) = 0.666666666666667⋅K⋅ρ                  
//   dρ                                             
//   
  dA_drho  = (2.0/3.0) * K * pow(t, -1.0/3.0);

  dfdt[0]  = dCs_dT * dT_dA * dA_drho;

# endif

  dfdt[0] *= sign;
  dfdt[0] *= LorentzFactor;

  return GSL_SUCCESS;
}


double Isentropic_Dens2Velocity ( double DensDown, struct Rarefaction *upstream )
{
  double VelyUp = upstream -> VelyUpStream;

  gsl_odeiv2_system sys = {func, jac, 1, &upstream};

  gsl_odeiv2_driver * d =  gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

  double t = 0.0, t1 = DensDown;

  double y[1] = { VelyUp };

  int NumStep = 100;

  for (int i = 1; i <= NumStep; i++)
  {
     double ti = i * t1 / (double)NumStep;

     int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

     if (status != GSL_SUCCESS)
     {
         printf ("error, return value=%d\n", status);
         break;
     }
  }

  gsl_odeiv2_driver_free (d);

  return y[0]; // return VelyDown;
}
