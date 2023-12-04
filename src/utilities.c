#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/struct.h"
#include "../include/prototypes.h"
#include "../include/global.h"
#include "../include/macro.h"


int GetAllInfomation( struct InitialCondition *IC, struct RiemannProblem *RP )
{
   struct Rarefaction Left;
   struct Rarefaction Right;

   double DensLeft      = IC -> DensLeft     ;
   double VelocityLeft  = IC -> VelocityLeft ;
   double PresLeft      = IC -> PresLeft     ;
   double DensRight     = IC -> DensRight    ;
   double VelocityRight = IC -> VelocityRight;
   double PresRight     = IC -> PresRight    ;

   int Pattern;

   Pattern = GetWavePattern( IC );

   double PresStar, VelocityStar;

   double up = 6e6;
   double lb = 6e5;


   PresStar = RootFinder( PresFunction, (void*)IC, 0.0, __DBL_EPSILON__, 0.5*(lb+up), lb, up, __FUNCTION__ );
   printf("PresStar=%20.16e\n", PresStar);

   double ShockVelocity_Left,  DensDown_Left;
   double ShockVelocity_Right, DensDown_Right;
   double HeadVelocity_Right,  TailVelocity_Right;
   double HeadVelocity_Left,   TailVelocity_Left;

   switch ( Pattern )
   {
      case 1:
         DensDown_Left = GetDensDown( PresLeft,  DensLeft,  PresStar );

         GetShockVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensDown_Left, &ShockVelocity_Left, NULL );

         DensDown_Right = GetDensDown( PresRight, DensRight, PresStar );

         GetShockVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, NULL, &ShockVelocity_Right );

         if ( ShockVelocity_Left > 0.0 )
            VelocityStar = GetVelocityDown( PresRight, DensRight, ShockVelocity_Right, PresStar, DensDown_Right );
         else
            VelocityStar = GetVelocityDown( PresLeft, DensLeft, ShockVelocity_Left, PresStar, DensDown_Left );

         RP->SS.Leftt.Right_Yes      = false;
         RP->SS.Leftt.ShockVelocity  = ShockVelocity_Left;
         RP->SS.Leftt.PresUpStream   = PresLeft;
         RP->SS.Leftt.DensUpStream   = DensLeft;
         RP->SS.Leftt.VelyUpStream   = VelocityLeft;
         RP->SS.Leftt.PresDownStream = PresStar;
         RP->SS.Leftt.DensDownStream = DensDown_Left;
         RP->SS.Leftt.VelyDownStream = VelocityStar;

         RP->SS.Right.Right_Yes      = true;
         RP->SS.Right.ShockVelocity  = ShockVelocity_Right;
         RP->SS.Right.PresUpStream   = PresRight;
         RP->SS.Right.DensUpStream   = DensRight;
         RP->SS.Right.VelyUpStream   = VelocityRight;
         RP->SS.Right.PresDownStream = PresStar;
         RP->SS.Right.DensDownStream = DensDown_Right;
         RP->SS.Right.VelyDownStream = VelocityStar;
         break;

      case 2:
         Left.Right_Yes      = false;
         Left.PresUpStream   = PresLeft;
         Left.DensUpStream   = DensLeft;
         Left.VelyUpStream   = VelocityLeft;
         Left.PresDownStream = PresStar;

         DensDown_Left = Isentropic_Pres2Dens( &Left );

         VelocityStar = Isentropic_Dens2Velocity( DensDown_Left, &Left );

         GetHeadTailVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensDown_Left, VelocityStar, &HeadVelocity_Left, &TailVelocity_Left, false );

         RP->RS.Leftt.Right_Yes      = false;
         RP->RS.Leftt.PresUpStream   = PresLeft;
         RP->RS.Leftt.DensUpStream   = DensLeft;
         RP->RS.Leftt.VelyUpStream   = VelocityLeft;
         RP->RS.Leftt.PresDownStream = PresStar;
         RP->RS.Leftt.DensDownStream = DensDown_Left;
         RP->RS.Leftt.VelyDownStream = VelocityStar;
         RP->RS.Leftt.VelocityHead   = HeadVelocity_Left;
         RP->RS.Leftt.VelocityTail   = TailVelocity_Left;

         DensDown_Right = GetDensDown( PresRight, DensRight, PresStar );

         GetShockVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, NULL, &ShockVelocity_Right );

         VelocityStar = GetVelocityDown( PresRight, DensRight, ShockVelocity_Right, PresStar, DensDown_Right );

         RP->RS.Right.Right_Yes      = true;
         RP->RS.Right.ShockVelocity  = ShockVelocity_Right;
         RP->RS.Right.PresUpStream   = PresRight;
         RP->RS.Right.DensUpStream   = DensRight;
         RP->RS.Right.VelyUpStream   = VelocityRight;
         RP->RS.Right.PresDownStream = PresStar;
         RP->RS.Right.DensDownStream = DensDown_Right;
         RP->RS.Right.VelyDownStream = VelocityStar;
         break;

      case 3:
         DensDown_Left = GetDensDown( PresLeft, DensLeft, PresStar );

         GetShockVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensDown_Left, &ShockVelocity_Left, NULL );

         VelocityStar = GetVelocityDown( PresLeft, DensLeft, ShockVelocity_Left, PresStar, DensDown_Left );

         RP->SR.Leftt.Right_Yes      = false;
         RP->SR.Leftt.ShockVelocity  = ShockVelocity_Left; // o
         RP->SR.Leftt.PresUpStream   = PresLeft; // o
         RP->SR.Leftt.DensUpStream   = DensLeft; // o
         RP->SR.Leftt.VelyUpStream   = VelocityLeft; // o
         RP->SR.Leftt.PresDownStream = PresStar; // o
         RP->SR.Leftt.DensDownStream = DensDown_Left; // o
         RP->SR.Leftt.VelyDownStream = VelocityStar; // o

         Right.Right_Yes      = true;
         Right.PresUpStream   = PresRight;
         Right.DensUpStream   = DensRight;
         Right.VelyUpStream   = VelocityRight;
         Right.PresDownStream = PresStar;

         DensDown_Right = Isentropic_Pres2Dens( &Right );

         VelocityStar = Isentropic_Dens2Velocity( DensDown_Right, &Right );

         GetHeadTailVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, VelocityStar, &HeadVelocity_Right, &TailVelocity_Right, true );

         RP->SR.Right.Right_Yes      = true;
         RP->SR.Right.PresUpStream   = PresRight;
         RP->SR.Right.DensUpStream   = DensRight;
         RP->SR.Right.VelyUpStream   = VelocityRight;
         RP->SR.Right.PresDownStream = PresStar;  // o
         RP->SR.Right.DensDownStream = DensDown_Right; // o
         RP->SR.Right.VelyDownStream = VelocityStar;       // x
         RP->SR.Right.VelocityHead   = HeadVelocity_Right; // x
         RP->SR.Right.VelocityTail   = TailVelocity_Right; // x
         break;

      case 4:
         Left.Right_Yes      = false;
         Left.PresUpStream   = PresLeft;
         Left.DensUpStream   = DensLeft;
         Left.VelyUpStream   = VelocityLeft;
         Left.PresDownStream = PresStar;

         DensDown_Left = Isentropic_Pres2Dens( &Left );

         VelocityStar = Isentropic_Dens2Velocity( DensDown_Left, &Left );

         GetHeadTailVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensDown_Left, VelocityStar, &HeadVelocity_Left, &TailVelocity_Left, false );

         RP->RR.Leftt.Right_Yes      = false;
         RP->RR.Leftt.PresUpStream   = PresLeft;
         RP->RR.Leftt.DensUpStream   = DensLeft;
         RP->RR.Leftt.VelyUpStream   = VelocityLeft;
         RP->RR.Leftt.PresDownStream = PresStar;
         RP->RR.Leftt.DensDownStream = DensDown_Left;
         RP->RR.Leftt.VelyDownStream = VelocityStar;
         RP->RR.Leftt.VelocityHead   = HeadVelocity_Left;
         RP->RR.Leftt.VelocityTail   = TailVelocity_Left;

         Right.Right_Yes      = true;
         Right.PresUpStream   = PresRight;
         Right.DensUpStream   = DensRight;
         Right.VelyUpStream   = VelocityRight;
         Right.PresDownStream = PresStar;

         DensDown_Right = Isentropic_Pres2Dens( &Right );

         VelocityStar = Isentropic_Dens2Velocity( DensDown_Right, &Right );

         GetHeadTailVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, VelocityStar, &HeadVelocity_Right, &TailVelocity_Right, true );

         RP->RR.Right.Right_Yes      = true;
         RP->RR.Right.PresUpStream   = PresRight;
         RP->RR.Right.DensUpStream   = DensRight;
         RP->RR.Right.VelyUpStream   = VelocityRight;
         RP->RR.Right.PresDownStream = PresStar;  // o
         RP->RR.Right.DensDownStream = DensDown_Right; // o
         RP->RR.Right.VelyDownStream = VelocityStar;       // x
         RP->RR.Right.VelocityHead   = HeadVelocity_Right; // x
         RP->RR.Right.VelocityTail   = TailVelocity_Right; // x
         break;

      default:
         printf("wave pattern was not found!!\n");
   } // switch ( Pattern )

   return Pattern;

} // FUNCTION : GetAllInfomation


double U2V( double U )
{
   return U / sqrt( 1.0 + U*U );
} // FUNCTION : U2V
