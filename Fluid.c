#include "Global.h"

double SoundSpeed( double Pres, double Dens )
{
	double Cs_2, Enthalpy;

	Enthalpy = Enthalpy( Pres, Dens );

	if( Pres > 0.0 )  
	{
	   Cs_2 = Gamma * Pres / Dens / Enthalpy;
	}
	else if ( Pres == 0.0 ) 
	{
	   Cs_2 = 0.0;
	}
	else
	{
	   printf("pressure == 0.0 was found!!\n");
	   exit(1);
	}

    return sqrt( Cs_2 );
}

double Enthalpy( double Pres, double Dens )
{
    return 1.0 + ( Gamma / Gamma_1 ) * ( Pres / Dens );
}

double TotalEngy ( double Pres, double Dens )
{
    double Enthalpy, Engy;

	Enthalpy = Enthalpy( Pres, Dens );
 
    Engy = Pres * ( Enthalpy / ( Enthalpy - 1.0 ) ) * ( Gamma / Gamma_1 ) - Pres;

    return Engy;
}



