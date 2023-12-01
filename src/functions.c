#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../includes/struct.h"
#include "../includes/prototypes.h"
#include "../includes/global.h"
#include "../includes/macro.h"

// Pattern 1: shock-shock
//         2: rarefaction-shock
//         3: shock-rarefaction
//         4: rarefaction-rarefaction


int GetWavePattern( struct InitialCondition *IC )
{
  double DensLeft      = IC -> DensLeft     ;
  double VelocityLeft  = IC -> VelocityLeft ;
  double PresLeft      = IC -> PresLeft     ;
  double DensRight     = IC -> DensRight    ;
  double VelocityRight = IC -> VelocityRight;
  double PresRight     = IC -> PresRight    ;

  double SS, RS, RR;
  double V_LC;
  bool Swap_Yes = false;

  bool Shock_Yes = true;


  // swap PresLeft and PresRight
  if ( PresLeft < PresRight )
  {
    PresLeft  = PresLeft * PresRight;
    PresRight = PresLeft / PresRight;
    PresLeft  = PresLeft / PresRight;

    DensLeft  = DensLeft * DensRight;
    DensRight = DensLeft / DensRight;
    DensLeft  = DensLeft / DensRight;

    double temp1, temp2;

	temp1 = -VelocityLeft;
	temp2 = -VelocityRight;

    VelocityLeft  = temp2;
    VelocityRight = temp1;

    Swap_Yes = true;
  }


  //===============================================
  // shock-shock
  double DensStarRight = GetDensDown( PresRight,  DensRight,  PresLeft );
  double V_RC = Velocity_RC (  PresLeft,  DensStarRight,  PresRight,  DensRight,  VelocityRight, Shock_Yes );
  V_LC = 0.0;
  RelativeVelocity( V_LC, V_RC, NULL, &SS );

  //===============================================
  // rarefaction-shock
  double PresStar = PresRight;

  struct Rarefaction RarefactionLeft;

  RarefactionLeft.PresUpStream   = PresLeft;
  RarefactionLeft.DensUpStream   = DensLeft;
  RarefactionLeft.VelyUpStream   = VelocityLeft;
  RarefactionLeft.PresDownStream = PresStar;
  RarefactionLeft.Right_Yes      = false;

  double DensStarLeft = Isentropic_Pres2Dens( &RarefactionLeft );
  double VelocityStar = Isentropic_Dens2Velocity( DensStarLeft, &RarefactionLeft );
  RelativeVelocity( VelocityLeft, VelocityStar, NULL, &V_LC );
  V_RC = 0.0;
  RelativeVelocity( V_LC, V_RC, NULL, &RS );

  RS = V_LC;

  //===============================================
  // rarefaction-rarefaction
  PresStar = 1e-5;

  RarefactionLeft.PresUpStream   = PresLeft;
  RarefactionLeft.DensUpStream   = DensLeft;
  RarefactionLeft.VelyUpStream   = VelocityLeft;
  RarefactionLeft.PresDownStream = PresStar;
  RarefactionLeft.Right_Yes      = false;

  DensStarLeft = Isentropic_Pres2Dens( &RarefactionLeft );
  VelocityStar = Isentropic_Dens2Velocity( DensStarLeft, &RarefactionLeft );
  RelativeVelocity( VelocityLeft, VelocityStar, NULL, &V_LC );


  struct Rarefaction RarefactionRight;

  RarefactionRight.PresUpStream   = PresRight;
  RarefactionRight.DensUpStream   = DensRight;
  RarefactionRight.VelyUpStream   = VelocityRight;
  RarefactionRight.PresDownStream = PresStar;
  RarefactionRight.Right_Yes      = true;

  DensStarRight = Isentropic_Pres2Dens( &RarefactionRight );
  VelocityStar = Isentropic_Dens2Velocity( DensStarRight, &RarefactionRight );
  RelativeVelocity( VelocityRight, VelocityStar, NULL, &V_RC );


  RelativeVelocity( V_LC, V_RC, NULL, &RR );

  //===============================================

  //relative velocity
  double VelocityLeftRight;
  int Pattern;

  RelativeVelocity( VelocityLeft, VelocityRight, NULL, &VelocityLeftRight );


  if ( VelocityLeftRight >= SS )
  {
    Pattern = 1;
	printf("SS pattern !!\n");
  }
  else if (  RS <= VelocityLeftRight && VelocityLeftRight < SS && Swap_Yes == false )
  {
    Pattern = 2;
	printf("RS pattern !!\n");
  }
  else if (  RS <= VelocityLeftRight && VelocityLeftRight < SS && Swap_Yes == true )
  {
    Pattern = 3;
	printf("SR pattern !!\n");
  }
  else if ( VelocityLeftRight < RS )
  {
    Pattern = 4;
	printf("RR pattern !!\n");
  }
  else
  {
    printf("wave pattern was not found!!\n");
	exit(1);
  }

  return Pattern;

}


double Velocity_LC ( double PresStar, double DensStarLeft, double PresLeft, double DensLeft, double VelocityLeft, bool Shock )
{
  double Velocity_LC;

  if ( Shock == true )
  {
     double EngyStarLeft, EngyLeft, EnthalpyStarLeft, TempStarLeft;

     EnthalpyStarLeft = GetEnthalpyDown( PresLeft, DensLeft, PresStar );

     TempStarLeft     = Enthalpy2Temperature( EnthalpyStarLeft );

     EngyStarLeft     = Flu_TotalInternalEngy ( PresStar, PresStar/TempStarLeft );

	 EngyLeft         = Flu_TotalInternalEngy( PresLeft, DensLeft );

	 // 4-velocity
     Velocity_LC  = ( PresStar     - PresLeft ) * ( EngyStarLeft - EngyLeft );
	 Velocity_LC /= ( EngyStarLeft + PresStar ) * ( EngyLeft     + PresLeft );

     if ( Velocity_LC < 0.0 )
	 {
	   printf("Velocity_LC = %20.16e !!\n", Velocity_LC);
	   exit(1);
	 }

	 Velocity_LC  = sqrt( Velocity_LC );
  }
  else
  {
     if ( DensStarLeft != DensStarLeft )
	 {
	   printf( "DensStarLeft should be provided!!\n" );
	   exit(1);
	 }


     double Velocity_C;
     struct Rarefaction upstream;

     upstream.Right_Yes    = false;
     upstream.DensUpStream = DensLeft;
     upstream.PresUpStream = PresLeft;
     upstream.VelyUpStream = VelocityLeft;


     Velocity_C  = Isentropic_Dens2Velocity( DensStarLeft, &upstream );
     RelativeVelocity( VelocityLeft, Velocity_C, NULL, &Velocity_LC );
  }

  return Velocity_LC;
}

double Velocity_RC ( double PresStar, double DensStarRight, double PresRight, double DensRight, double VelocityRight, bool Shock )
{
  double Velocity_RC;

  if ( Shock == true )
  {
     double EngyStarRight, EngyRight, EnthalpyStarRight, TempStarRight;

     EnthalpyStarRight = GetEnthalpyDown( PresRight, DensRight, PresStar );

     TempStarRight     = Enthalpy2Temperature( EnthalpyStarRight );

     EngyStarRight     = Flu_TotalInternalEngy( PresStar, PresStar/TempStarRight );

	 EngyRight         = Flu_TotalInternalEngy( PresRight, DensRight );

	 // 4-velocity
     Velocity_RC  = ( PresStar      - PresRight ) * ( EngyStarRight - EngyRight );
	 Velocity_RC /= ( EngyStarRight + PresStar  ) * ( EngyRight     + PresRight );

     if ( Velocity_RC < 0.0 )
	 {
	   printf("Velocity_RC = %20.16e\n !!", Velocity_RC);
	   exit(1);
	 }

	 Velocity_RC  = -sqrt( Velocity_RC );
  }
  else
  {
     if ( DensStarRight != DensStarRight )
	 {
	   printf( "DensStarRight should be provided!!\n" );
	   exit(1);
	 }


     double Velocity_C;
     struct Rarefaction upstream;

     upstream.Right_Yes    = true;
     upstream.DensUpStream = DensRight;
     upstream.PresUpStream = PresRight;
     upstream.VelyUpStream = VelocityRight;


     Velocity_C  = Isentropic_Dens2Velocity( DensStarRight, &upstream );
     RelativeVelocity( VelocityRight, Velocity_C, NULL, &Velocity_RC );
  }

  return Velocity_RC;
}


double PresFunction( double PresStar, void  *params )
{

  if ( PresStar <= 0.0 )
  {
    printf("PresStar=%20.16e !!\n", PresStar);
	exit(1);
  }

  struct InitialCondition *IC = ( struct InitialCondition * ) params;
  double DensLeft      = IC -> DensLeft     ;
  double VelocityLeft  = IC -> VelocityLeft ;
  double PresLeft      = IC -> PresLeft     ;
  double DensRight     = IC -> DensRight    ;
  double VelocityRight = IC -> VelocityRight;
  double PresRight     = IC -> PresRight    ;

  bool Shock_Yes = true;
  bool Shock_No  = false;
  double V_LC, V_RC, V_LR;
  double DensStarLeft, DensStarRight;


  if ( PresStar >= MAX(PresLeft, PresRight) )
  {
    V_LC = Velocity_LC( PresStar, NAN, PresLeft,   DensLeft, VelocityLeft, Shock_Yes ); // left side of eq. (4.161)
    V_RC = Velocity_RC( PresStar, NAN, PresRight, DensRight, VelocityRight, Shock_Yes ); // right side of eq. (4.161)

	V_LR = - V_RC * sqrt(1.0 + V_LC * V_LC) + sqrt(1.0 + V_RC * V_RC) * V_LC;
  }
  else if ( MIN(PresLeft, PresRight) <= PresStar && PresStar < MAX(PresLeft, PresRight) && PresLeft >= PresRight  )
  {
    struct Rarefaction rarefation;
    rarefation.DensUpStream   = DensLeft;
    rarefation.PresUpStream   = PresLeft;
    rarefation.PresDownStream = PresStar;

    DensStarLeft = Isentropic_Pres2Dens( &rarefation );

    V_LC = Velocity_LC( PresStar, DensStarLeft, PresLeft,   DensLeft, VelocityLeft, Shock_No  ); // eq. (4.168)
    V_RC = Velocity_RC( PresStar, NAN,         PresRight,  DensRight, VelocityRight,  Shock_Yes ); // right side of eq. (4.161)

	V_LR = - V_RC * sqrt(1.0 + V_LC * V_LC) + sqrt(1.0 + V_RC * V_RC) * V_LC;
  }
  else if ( MIN(PresLeft, PresRight) <= PresStar && PresStar < MAX(PresLeft, PresRight) && PresLeft <= PresRight )
  {
    struct Rarefaction rarefation;
    rarefation.DensUpStream   = DensRight;
    rarefation.PresUpStream   = PresRight;
    rarefation.PresDownStream = PresStar;

    DensStarRight = Isentropic_Pres2Dens( &rarefation );

    V_LC = Velocity_LC( PresStar, NAN,          PresLeft,  DensLeft, VelocityLeft, Shock_Yes );
    V_RC = Velocity_RC( PresStar, DensStarRight, PresRight, DensRight, VelocityRight,  Shock_No  );

	V_LR = - V_RC * sqrt(1.0 + V_LC * V_LC) + sqrt(1.0 + V_RC * V_RC) * V_LC;
  }
  else if ( PresStar < MIN(PresLeft, PresRight)  )
  {
    struct Rarefaction Left;
    struct Rarefaction Right;

    Left.Right_Yes       = false;
    Left.PresUpStream    = PresLeft;
    Left.DensUpStream    = DensLeft;
    Left.VelyUpStream    = VelocityLeft;
    Left.PresDownStream  = PresStar;

    Right.Right_Yes      = true;
    Right.PresUpStream   = PresRight;
    Right.DensUpStream   = DensRight;
    Right.VelyUpStream   = VelocityRight;
    Right.PresDownStream = PresStar;

    DensStarLeft         = Isentropic_Pres2Dens( &Left );
    DensStarRight        = Isentropic_Pres2Dens( &Right );

    V_LC = Velocity_LC( PresStar, DensStarLeft, PresLeft,   DensLeft, VelocityLeft, Shock_No  ); // eq. (4.168)
    V_RC = Velocity_RC( PresStar, DensStarRight, PresRight, DensRight, VelocityRight,  Shock_No  );
    RelativeVelocity( V_LC, V_RC, NULL, &V_LR );
  }else REPORT_ERROR;

  double VelocityLeftRight;


  VelocityLeftRight = - VelocityRight * sqrt(1.0 + VelocityLeft *VelocityLeft )
	                  +  VelocityLeft * sqrt(1.0 + VelocityRight*VelocityRight);

  return VelocityLeftRight - V_LR;
}



void QuadraticSolver( double A, double B, double C , double *PlusRoot, double *MinusRoot)
{
  double Delta;

  Delta = sqrt( B*B - 4.0*A*C );

  if ( Delta != Delta )
  {
    printf("no solutions in quadratic equation!!\n");
	exit(1);
  }

  if ( PlusRoot  != NULL )  *PlusRoot  = -2.0*C/( +B + Delta);
  if ( MinusRoot != NULL )  *MinusRoot = +2.0*C/( -B + Delta);

}

// U_ab is the relative 4-velocity of U_a w.r.t. Ub

void RelativeVelocity( double Ua, double Ub, double *LorentzFactor_ab, double *U_ab )
{
  double LorentzFactor_a = sqrt(1.0 + Ua*Ua);
  double LorentzFactor_b = sqrt(1.0 + Ub*Ub);


  if ( U_ab != NULL )
  *U_ab             =  -Ub*LorentzFactor_a + LorentzFactor_b*Ua;

  if ( LorentzFactor_ab != NULL )
  *LorentzFactor_ab =  LorentzFactor_a*LorentzFactor_b - Ua*Ub;
}
