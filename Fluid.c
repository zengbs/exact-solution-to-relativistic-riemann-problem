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
       CsSqr *= 5.0*h - 8.0*T;
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
    return 2.5*Temp + sqrt( 2.25*Temp*Temp + 1.0 )
#   endif
}


double Flu_TotalEngy ( double Pres, double Dens )
{
    double Enthalpy, Engy;
#   if ( EOS == GAMMA )
	Enthalpy = Flu_Enthalpy( Pres, Dens );
 
    Engy = Pres * ( Enthalpy / ( Enthalpy - 1.0 ) ) * ( Gamma / Gamma_1 ) - Pres;
#   endif
    return Engy;
}



