#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"
#include "Macro.h"


// 4-sound speed
// TM EoS:
//          ⎛               _____________⎞
//          ⎜              ╱       2     ⎟
//        T⋅⎝4.5⋅T + 5.0⋅╲╱  2.25⋅T  + 1 ⎠
//   ──────────────────────────────────────
//                      _____________
//         2           ╱       2
//   13.5⋅T  + 7.0⋅T⋅╲╱  2.25⋅T  + 1  + 3.0


double Flu_SoundSpeed( double Temp )
{
   double CsSqr;

   if ( Temp > 0.0 )
   {
#     if ( EOS == GAMMA )
      CsSqr  = Gamma * Temp * Gamma_1;
      CsSqr /= Gamma * Temp * (2.0 - Gamma) + Gamma_1;
#     elif ( EOS == TM )
      CsSqr  = Temp * ( 4.5*Temp + 5.0*sqrt( 2.25*Temp*Temp + 1.0 ) );
      CsSqr /= 13.5*Temp*Temp + 7.0*Temp*sqrt( 2.25*Temp*Temp + 1.0 ) + 3.0;
#     endif
   } else
   {
      printf("Temp = %e !!\n", Temp);
      exit(1);
   }

   return sqrt( CsSqr );
} // FUNCTION : Flu_SoundSpeed



double Flu_Enthalpy( double Pres, double Dens )
{
#  if ( EOS == GAMMA )
   return 1.0 + ( Gamma / Gamma_1 ) * ( Pres / Dens );
#  elif ( EOS == TM )
   double Temp = Pres/Dens;
   return 2.5*Temp + sqrt( 2.25*Temp*Temp + 1.0 );
#  endif
} // FUNCTION : Flu_Enthalpy



// h = e/rho + p/rho
// e = rho*h - p
double Flu_TotalInternalEngy( double Pres, double Dens )
{
   double Engy;
   double Temp = Pres/Dens;

#  if ( EOS == GAMMA )
   Engy = Dens * ( 1.0 + Temp/Gamma_1 );
#  elif ( EOS == TM )
   Engy = Dens * ( 1.5*Temp + sqrt( 2.25*Temp*Temp + 1.0 ) );
#  endif

   return Engy;
} // FUNCTION : Flu_TotalInternalEngy



double Enthalpy2Temperature( double Enthalpy )
{
   double Temp;
   double H_Tilde = Enthalpy - 1.0;
#  if ( EOS == GAMMA )
   Temp  = Gamma_1/Gamma;
   Temp *= H_Tilde;
#  elif ( EOS == TM )
   Temp  = 2.0*H_Tilde*H_Tilde + 4.0*H_Tilde;
   Temp /= 5.0*(H_Tilde+1.0)   + sqrt( 9.0*H_Tilde*H_Tilde+18.0*H_Tilde+25.0 );
#  endif
   return Temp;
} // FUNCTION : Enthalpy2Temperature
