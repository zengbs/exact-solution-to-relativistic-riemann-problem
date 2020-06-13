#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"
#include "Macro.h"

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
  double V_LC, V_RC;
  bool Shock_Yes = true;
  bool Shock_No  = false;
  bool Swap_Yes = false;

  double A_PlusLeft, A_MinusLeft, A_PlusRight, A_MinusRight;

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
  double A, B, C, EnthalpyRight, Root, Engy_Temp, EngyRight;

  EnthalpyRight = Flu_Enthalpy ( PresRight, DensRight );
  EngyRight     = Flu_TotalInternalEngy( PresRight, DensRight );

  A = 1.0 + ( Gamma_1 / Gamma ) * ( PresRight / PresLeft - 1.0 );
  B = - ( Gamma_1 / Gamma ) * ( PresRight / PresLeft - 1.0 );
  C = EnthalpyRight * ( PresRight - PresLeft ) / DensRight - SQR(EnthalpyRight);


  QuadraticSolver( A, B, C, &Root, NULL ); 

  Engy_Temp = ( Gamma / Gamma_1 ) * ( Root / (Root-1.0) ) * PresLeft - PresLeft;// eq. (4.165)

  // 4-velocity
  SS = sqrt( ( Engy_Temp - EngyRight )*( PresLeft - PresRight )/( EngyRight + PresRight )/( Engy_Temp + PresLeft )  );

  //===============================================
  // rarefaction-shock
  double DensStarLeft = DensLeft*pow(  PresRight/PresLeft, 1.0/Gamma );

  A_PlusRight = A_PlusFun( PresRight / DensStarLeft );

  A_PlusLeft  = A_PlusFun( PresLeft  / DensLeft    );

  V_LC  = ( A_PlusRight - A_PlusLeft )/sqrt( 4.0*A_PlusRight*A_PlusLeft );

  // 4-velocity
  RS = V_LC;

  //===============================================
  // rarefaction-rarefaction

  A_PlusLeft  = A_PlusFun ( PresLeft  / DensLeft  );
  A_MinusRight = A_MinusFun( PresRight / DensRight );

  // 4-velocity
  RR = ( A_MinusRight - A_PlusLeft )/sqrt( 4.0 * A_PlusLeft * A_MinusRight );

  // relative velocity  
  double RelitiveVelocity;
  int Pattern;


  RelitiveVelocity = - VelocityRight * sqrt(1.0 + VelocityLeft *VelocityLeft )
		             +  VelocityLeft * sqrt(1.0 + VelocityRight*VelocityRight);


  if ( RelitiveVelocity >= SS )
  {
    Pattern = 1;
	//printf("you have shock-shock wave pattern !!\n");
  }
  else if (  RS <= RelitiveVelocity && RelitiveVelocity < SS && Swap_Yes == false )
  {
    Pattern = 2;
	//printf("you have rarefaction-shock wave pattern !!\n");
  }
  else if (  RS <= RelitiveVelocity && RelitiveVelocity < SS && Swap_Yes == true )
  {
    Pattern = 3;
	//printf("you have shock-rarefaction wave pattern !!\n");
  }
  else if ( RR <= RelitiveVelocity && RelitiveVelocity < RS )
  {
    Pattern = 4;
	//printf("you have rarefaction-rarefaction wave pattern !!\n");
  }
  else
  {
    printf("wave pattern was not found!!\n");
	exit(1);
  }

  return Pattern;

}


double Velocity_LC ( double PresStar, double DensStarLeft, double PresLeft, double DensLeft, bool Shock )
{
  double Velocity_LC;

  if ( Shock == true )
  {
     double EngyStarLeft, EngyLeft, EnthalpyStarLeft, TempStarLeft;

     EnthalpyStarLeft = TaubAdiabatic( PresLeft, DensLeft, PresStar );

     TempStarLeft     = Enthalpy2Temperature( EnthalpyStarLeft );

     EngyStarLeft     = Flu_TotalInternalEngy ( PresStar, PresStar/TempStarLeft );

	 EngyLeft         = Flu_TotalInternalEngy( PresLeft, DensLeft );

	 // 4-velocity
     Velocity_LC  = ( PresStar     - PresLeft ) * ( EngyStarLeft - EngyLeft );
	 Velocity_LC /= ( EngyStarLeft + PresStar ) * ( EngyLeft     + PresLeft );



     if ( Velocity_LC < 0.0 )
	 {
	   printf("Velocity_LC = %e\n !!", Velocity_LC);
	   exit(1);
	 }

	 Velocity_LC  = sqrt( Velocity_LC );

  
     return Velocity_LC;
  }
  else
  {
     if ( DensStarLeft != DensStarLeft )
	 {
	   printf( "DensStarLeft should be provided!!\n" );
	   exit(1);
	 }


	 // 4-velocity
	 double A_PlusStar, A_PlusLeft;

	 A_PlusStar  = A_PlusFun( PresStar / DensStarLeft);
	 A_PlusLeft  = A_PlusFun( PresLeft / DensLeft );

     Velocity_LC = ( A_PlusStar - A_PlusLeft )/sqrt( 4.0 * A_PlusStar * A_PlusLeft );

     return Velocity_LC; 
  }
}

double Velocity_RC ( double PresStar, double DensStarRight, double PresRight, double DensRight, bool Shock )
{
  double Velocity_RC;

  if ( Shock == true )
  {
     double EngyStarRight, EngyRight, EnthalpyStarRight, TempStarRight;

     EnthalpyStarRight = TaubAdiabatic( PresRight, DensRight, PresStar );

     TempStarRight     = Enthalpy2Temperature( EngyStarRight );

     EngyStarRight     = Flu_TotalInternalEngy( PresStar, PresStar/TempStarRight );

	 EngyRight         = Flu_TotalInternalEngy( PresRight, DensRight );


	 // 4-velocity
     Velocity_RC  = ( PresStar      - PresRight ) * ( EngyStarRight - EngyRight );
	 Velocity_RC /= ( EngyStarRight + PresStar  ) * ( EngyRight     + PresRight );

     if ( Velocity_RC < 0.0 )
	 {
	   printf("Velocity_RC = %e\n !!", Velocity_RC);
	   exit(1);
	 }

	 Velocity_RC  = -sqrt( Velocity_RC );

     return Velocity_RC;
  }
  else
  {
     if ( DensStarRight != DensStarRight )
	 {
	   printf( "DensStarRight should be provided!!\n" );
	   exit(1);
	 }


	 // 4-velocity
	 double A_MinusStar, A_MinusRight;

	 A_MinusStar  = A_MinusFun( PresStar / DensStarRight);
	 A_MinusRight  = A_MinusFun( PresRight / DensRight );

     Velocity_RC = ( A_MinusStar - A_MinusRight )/sqrt( 4.0 * A_MinusStar * A_MinusRight );
  
     return Velocity_RC; 
  }
}

double A_PlusFun ( double Temp )
{
    double Sqrt_Gamma_1 = sqrt( Gamma_1 );
    double A_Plus;

    A_Plus  = SQR( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );
	A_Plus /= Gamma_1;
	A_Plus  = pow(A_Plus, +2.0/Sqrt_Gamma_1);

	return A_Plus;
}

double A_MinusFun ( double Temp )
{
    double Sqrt_Gamma_1 = sqrt( Gamma_1 );
    double A_Minus;

    A_Minus  = SQR( sqrt( Gamma_1 + Gamma * Temp ) + sqrt( Gamma * Temp ) );
	A_Minus /= Gamma_1;
	A_Minus  = pow(A_Minus, -2.0/Sqrt_Gamma_1);

	return A_Minus;
}


double PresFunction( double PresStar, void  *params )
{

  if ( PresStar <= 0.0 )
  {
    printf("PresStar=%e !!\n", PresStar);
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
    V_LC = Velocity_LC( PresStar, NAN, PresLeft,   DensLeft, Shock_Yes ); // left side of eq. (4.161)
    V_RC = Velocity_RC( PresStar, NAN, PresRight, DensRight, Shock_Yes ); // right side of eq. (4.161)

	V_LR = - V_RC * sqrt(1.0 + V_LC * V_LC) + sqrt(1.0 + V_RC * V_RC) * V_LC;
  }
  else if ( MIN(PresLeft, PresRight) <= PresStar && PresStar < MAX(PresLeft, PresRight) && PresLeft >= PresRight  )
  {
    DensStarLeft = DensLeft*pow(  PresStar/PresLeft, 1.0/Gamma );

    V_LC = Velocity_LC( PresStar, DensStarLeft, PresLeft,   DensLeft, Shock_No  ); // eq. (4.168)
    V_RC = Velocity_RC( PresStar, NAN,         PresRight,  DensRight, Shock_Yes ); // right side of eq. (4.161) 

	V_LR = - V_RC * sqrt(1.0 + V_LC * V_LC) + sqrt(1.0 + V_RC * V_RC) * V_LC;
  }
  else if ( MIN(PresLeft, PresRight) <= PresStar && PresStar < MAX(PresLeft, PresRight) && PresLeft <= PresRight )
  {
    DensStarRight = DensRight*pow(  PresStar/PresRight, 1.0/Gamma );
  
    V_LC = Velocity_LC( PresStar, NAN,          PresLeft,  DensLeft, Shock_Yes );
    V_RC = Velocity_RC( PresStar, DensStarRight, PresRight, DensRight, Shock_No  );

	V_LR = - V_RC * sqrt(1.0 + V_LC * V_LC) + sqrt(1.0 + V_RC * V_RC) * V_LC;
  }
  else if ( PresStar < MIN(PresLeft, PresRight)  )
  {
    DensStarLeft  = DensLeft *pow(  PresStar/PresLeft,  1.0/Gamma );
    DensStarRight = DensRight*pow(  PresStar/PresRight, 1.0/Gamma );

	double V_LC, V_RC;
    double A_PlusStar, A_MinusStar;
    double A_PlusLeft, A_MinusRight;

    A_PlusStar   = A_PlusFun( PresStar  / DensStarLeft );
    A_PlusLeft   = A_PlusFun( PresLeft  / DensLeft     );

    A_MinusStar  = A_MinusFun( PresStar  / DensStarRight );
    A_MinusRight = A_MinusFun( PresRight / DensRight    );

    V_LC  = A_PlusStar - A_PlusLeft;
    V_LC /= sqrt(4.0 * A_PlusStar * A_PlusLeft);

    V_RC  = A_MinusStar - A_MinusRight;
    V_RC /= sqrt(4.0 * A_MinusStar * A_MinusRight);

	V_LR = -sqrt(1.0+V_LC*V_LC)*V_RC + sqrt(1.0+V_RC*V_RC)*V_LC;
  }

  double RelitiveVelocity;


  RelitiveVelocity = - VelocityRight * sqrt(1.0 + VelocityLeft *VelocityLeft )
		             +  VelocityLeft * sqrt(1.0 + VelocityRight*VelocityRight);

  return RelitiveVelocity - V_LR;
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
