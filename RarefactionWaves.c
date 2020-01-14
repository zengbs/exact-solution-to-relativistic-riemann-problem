#include <stdbool.h>
#include <math.h>
#include "Global.h"
#include "Prototypes.h"

void GetHeadTailVelocity( double PresHead, double DensHead, double VelocityHead,
			              double PresTail, double DensTail, double VelocityTail,
                          double *HeadVelocity, double *TailVelocity, bool Right_Yes )
{
  double Cs_Head, Cs_Tail;

  Cs_Head = Flu_SoundSpeed( PresHead, DensHead ); 
  Cs_Tail = Flu_SoundSpeed( PresTail, DensTail ); 

  if ( Right_Yes )
  {
     *HeadVelocity = ( VelocityHead + Cs_Head ) / ( 1.0 + VelocityHead*Cs_Head );
     *TailVelocity = ( VelocityTail + Cs_Tail ) / ( 1.0 + VelocityTail*Cs_Tail );
  }
  else
  {
     *HeadVelocity = ( VelocityHead - Cs_Head ) / ( 1.0 - VelocityHead*Cs_Head );
     *TailVelocity = ( VelocityTail - Cs_Tail ) / ( 1.0 - VelocityTail*Cs_Tail );
  }
}


void GetDensInFan( double Cs2, double PresHead, double DensHead, double Xi )
{
  double k = PresHead * pow (DensHead, -Gamma);
  double tmp = ( 1.0 / (Cs2 * Cs2)) - (1.0 / Gamma_1 );
  double Dens = pow ( k * Gamma * tmp, -1.0 / Gamma_1 );

  return Dens;
}

void GetPresInFan( double DensInFan, double PresHead, double DensHead )
{
  double Pres = PresHead * pow (DensInFan / DensHead, Gamma);

  return Pres;
}

void GetVelicityInFan( double Cs2, double Xi )
{
  double Velocity = (Cs2 + Xi) / (1.0 + Cs2 * Xi);

  return Velocity;
}


double SoundSpeedFunction ( double Cs, void *params )
{
  struct RareFaction *Fan = ( struct RareFaction * ) params;

  bool   Right_Yes    = Fan -> Right_Yes   ;
  double PresHead     = Fan -> PresHead    ;
  double DensHead     = Fan -> DensHead    ;
  double VelocityHead = Fan -> VelocityHead;
  double PresTail     = Fan -> PresTail    ;
  double DensTail     = Fan -> DensTail    ;
  double VelocityTail = Fan -> VelocityTail;
  double Xi           = Fan -> Xi          ;

 
  double Velocity, Var0, Var1, Cs_Head;

  double Sqrt_Gamma_1 = sqrt(Gamma_1);

  Cs_Head = Flu_SoundSpeed( PresHead, DensHead );

  if ( Right_Yes )
  {
    Velocity = ( Xi - Cs )/( 1.0 - Cs * Xi );

    Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	Var0  = pow( Var0, -2.0 / Sqrt_Gamma_1 );

	Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity );


	Var1  = ( Sqrt_Gamma_1 + Cs_Head ) / ( Sqrt_Gamma_1 - Cs_Head );

	Var1  = pow( Var1, -2.0 / Sqrt_Gamma_1 );

	Var1 *= ( 1.0 + VelocityHead ) / ( 1.0 - VelocityHead );
  }
  else
  {
    Velocity = ( Xi + Cs )/( 1.0 + Cs * Xi );

    Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	Var0  = pow( Var0, +2.0 / Sqrt_Gamma_1 );

	Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity );


	Var1  = ( Sqrt_Gamma_1 + Cs_Head ) / ( Sqrt_Gamma_1 - Cs_Head );

	Var1  = pow( Var1, -2.0 / Sqrt_Gamma_1 );

	Var1 *= ( 1.0 + VelocityHead ) / ( 1.0 - VelocityHead );
  }

  return Var1 - Var0;

}

double GetSoundSpeedInFan ( struct RareFaction *Fan )
{
  double Cs2;

  Cs2 = RootFinder( SoundSpeedFunction, (void*)Fan, 0.0, __DBL_EPSILON__, 0.5, 1e-2, 1.0 );

  return Cs2;
}

