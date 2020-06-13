#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"


double Flu_SoundSpeed( double Temp )
{
	double CsSqr;

	if( Temp > 0.0 )  
	{
#      if ( EOS == GAMMA )
	   CsSqr  = Gamma * Temp * Gamma_1;
	   CsSqr /= Gamma * Temp * (2.0 - Gamma) + Gamma_1;
#      elif ( EOS == TM )
       double h = Flu_Enthalpy(Temp, 1.0);
       CsSqr  = Temp / (3.0*h);
       CsSqr *= 5.0*h - 8.0*Temp;
       CsSqr /= h - Temp;
#      endif
	}
	else
	{
	   printf("Temp = %e !!\n", Temp);
	   exit(1);
	}

    return sqrt( CsSqr );
}


double Flu_Enthalpy( double Pres, double Dens )
{
#   if ( EOS == GAMMA )
    return 1.0 + ( Gamma / Gamma_1 ) * ( Pres / Dens );
#   elif ( EOS == TM )
    double Temp = Pres/Dens;
    return 2.5*Temp + sqrt( 2.25*Temp*Temp + 1.0 );
#   endif
}



// h = e/rho + p/rho
// e = rho*h - p
double Flu_TotalInternalEngy ( double Pres, double Dens )
{
    double Enthalpy, Engy;

	Enthalpy = Flu_Enthalpy( Pres, Dens );
    Engy = Dens*Enthalpy - Pres; 

    return Engy;
}


double Enthalpy2Temperature( double Enthalpy )
{
  double Temp;
  double H_Tilde = Enthalpy - 1.0;
# if ( EOS == GAMMA )
  Temp  = (Gamma-1.0)/Gamma;
  Temp *= H_Tilde;
# elif ( EOS == TM )
  Temp  = 2.0*H_Tilde*H_Tilde + 4.0*H_Tilde;
  Temp /= 5.0*(H_Tilde+1.0) + sqrt( 9.0*H_Tilde*H_Tilde+18.0*H_Tilde+25.0 );
# endif
  return Temp;
}
