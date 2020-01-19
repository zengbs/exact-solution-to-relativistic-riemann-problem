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
	   CsSqr  = Gamma * Temp * Gamma_1;
	   CsSqr /= Gamma * Temp * (2.0 - Gamma) + Gamma_1;
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
    return 1.0 + ( Gamma / Gamma_1 ) * ( Pres / Dens );
}

double Flu_TotalEngy ( double Pres, double Dens )
{
    double Enthalpy, Engy;

	Enthalpy = Flu_Enthalpy( Pres, Dens );
 
    Engy = Pres * ( Enthalpy / ( Enthalpy - 1.0 ) ) * ( Gamma / Gamma_1 ) - Pres;

    return Engy;
}



