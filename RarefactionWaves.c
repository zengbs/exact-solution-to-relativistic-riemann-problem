void FanVelocity( double PresHead, double DensHead, double VelocityHead,
				  double PresTail, double DensTail, double VelocityTail,
                  double *HeadVelocity, double *TailVelocity, bool Right_Yes )
{
  double Cs_Head, Cs_Tail;

  Cs_Head = Flu_SoundSpeeed( PresHead, DensHead ); 
  Cs_Tail = Flu_SoundSpeeed( PresTail, DensTail ); 

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




void RareFactionFan( double PresHead, double DensHead, double VelocityHead,
                     double PresTail, double DensTail, double VelocityTail, )

{
                   




}

void GetSoundSpeedInFan ( double Cs, void *params )
{
  struct RareFactionFan *Fan = ( struct RareFactionFan * ) params;

  bool   Right_Yes    = Fan -> Right_Yes   ;
  double PresHead     = Fan -> PresHead    ;
  double DensHead     = Fan -> DensHead    ;
  double VelocityHead = Fan -> VelocityHead;
  double PresTail     = Fan -> PresTail    ;
  double DensTail     = Fan -> DensTail    ;
  double VelocityTail = Fan -> VelocityTail;
  double Xi           = Fan -> Xi          ;

 
  double Velocity, Var0, Var1, Cs_Head;

  Cs_Head = Flu_SoundSpeed( PresHead, DensHead );

  if ( Right_Yes )
  {
    Velocity = ( Xi - Cs )/( 1.0 - Cs * Xi );

    Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	Var0  = pow( Var0, -2.0 / Sqrt_Gamma_1 );

	Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity )


	Var1  = ( Sqrt_Gamma_1 + Cs_Head ) / ( Sqrt_Gamma_1 - Cs_Head );

	Var1  = pow( Var1, -2.0 / Sqrt_Gamma_1 );

	Var1 *= ( 1.0 + VelocityHead ) / ( 1.0 - VelocityHead );
  }
  else
  {
    Velocity = ( Xi + Cs )/( 1.0 + Cs * Xi );

    Var0  = ( Sqrt_Gamma_1 + Cs ) / ( Sqrt_Gamma_1 - Cs );

	Var0  = pow( Var0, +2.0 / Sqrt_Gamma_1 );

	Var0 *= ( 1.0 + Velocity ) / ( 1.0 - Velocity )


	Var1  = ( Sqrt_Gamma_1 + Cs_Head ) / ( Sqrt_Gamma_1 - Cs_Head );

	Var1  = pow( Var1, -2.0 / Sqrt_Gamma_1 );

	Var1 *= ( 1.0 + VelocityHead ) / ( 1.0 - VelocityHead );
  }

  return Var1 - Var0;

}
