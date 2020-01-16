#include <stdbool.h>
#include <stdio.h>
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

double GetDensDownRarefaction( double PresDown, double PresUp, double DensUp )
{
  double DensDown;

  DensDown = pow( DensUp, Gamma )*PresDown / PresUp;

  DensDown = pow( DensDown, 1.0 / Gamma );
 
  return DensDown;
}


double GetVelocityDownRarefaction( double PresDown, double DensDown, double PresUp, double DensUp, double VelocityUp )
{
  double Velocity;

  Velocity  = (1.0 + VelocityUp) * A_PlusFun( PresDown, DensDown, PresUp, DensUp ) - (1.0 - VelocityUp);
  Velocity /= (1.0 + VelocityUp) * A_PlusFun( PresDown, DensDown, PresUp, DensUp ) + (1.0 - VelocityUp);

  return Velocity;
}

double GetDensInFan( double Cs, double PresUp, double DensUp )
{
  double k = PresUp * pow (DensUp, -Gamma);
  double tmp = ( 1.0 / (Cs * Cs)) - (1.0 / Gamma_1 );
  double Dens = pow ( k * Gamma * tmp, -1.0 / Gamma_1 );

  return Dens;
}

double GetPresInFan( double DensInFan, double PresUp, double DensUp )
{
  double Pres = PresUp * pow (DensInFan / DensUp, Gamma);

  return Pres;
}

double GetVelocityInFan( double Cs2, double Xi )
{
  double Velocity = (Cs2 + Xi) / (1.0 + Cs2 * Xi);

  return Velocity;
}


double GetSoundSpeedInFan ( struct Rarefaction *Rarefaction )
{
  double Cs;

  Cs = RootFinder( SoundSpeedFunction, (void*)Rarefaction, 0.0, __DBL_EPSILON__, 0.11, 0.0, sqrt (Gamma_1) - 1e-16 );

  return Cs;
}

double SoundSpeedFunction ( double Cs, void *params )
{
  struct Rarefaction *Fan = ( struct Rarefaction * ) params;

  bool   Right_Yes    = Fan -> Right_Yes   ;
  double PresUp       = Fan -> PresUpStream;
  double DensUp       = Fan -> DensUpStream;
  double VelocityUp   = Fan -> VelyUpStream;
  double Xi           = Fan -> Xi          ;
 
  double Velocity, Var0, Var1, Cs_Up;

  double Sqrt_Gamma_1 = sqrt(Gamma_1);

  Cs_Up = Flu_SoundSpeed( PresUp, DensUp );

  if ( Right_Yes )
  {
    Velocity = ( Xi - Cs )/( 1.0 - Cs * Xi );

    Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	Var0  = pow( Var0, -2.0 / Sqrt_Gamma_1 );

	Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity );


	Var1  = ( Sqrt_Gamma_1 + Cs_Up ) / ( Sqrt_Gamma_1 - Cs_Up );

	Var1  = pow( Var1, -2.0 / Sqrt_Gamma_1 );

	Var1 *= ( 1.0 + VelocityUp ) / ( 1.0 - VelocityUp );
  }
  else
  {
    Velocity = ( Xi + Cs )/( 1.0 + Cs * Xi );

    Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	Var0  = pow( Var0, +2.0 / Sqrt_Gamma_1 );

	Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity );


	Var1  = ( Sqrt_Gamma_1 + Cs_Up ) / ( Sqrt_Gamma_1 - Cs_Up );

	Var1  = pow( Var1, +2.0 / Sqrt_Gamma_1 );

	Var1 *= ( 1.0 + VelocityUp ) / ( 1.0 - VelocityUp );
  }
  return Var1 - Var0;

}

