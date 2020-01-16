#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Struct.h"
#include "Prototypes.h"
#include "Global.h"
#include "Macro.h"


int GetAllInfomation( struct InitialCondition *IC, struct RiemannProblem *RP )
{
  double DensLeft      = IC -> DensLeft     ;
  double VelocityLeft  = IC -> VelocityLeft ;
  double PresLeft      = IC -> PresLeft     ;
  double DensRight     = IC -> DensRight    ;
  double VelocityRight = IC -> VelocityRight;
  double PresRight     = IC -> PresRight    ;

  int Pattern;

  Pattern = GetWavePattern( IC );
  
  double PresStar, VelocityStar;

  PresStar = RootFinder( PresFunction, (void*)IC, 0.0, __DBL_EPSILON__, 0.5, 1e-2, 50.0 );


  double ShockVelocity_Left,  DensDown_Left;
  double ShockVelocity_Right, DensDown_Right;
  double HeadVelocity_Right,  TailVelocity_Right;
  double HeadVelocity_Left,   TailVelocity_Left;

  switch ( Pattern )
  {
	 case 1:
         DensDown_Left  = GetDensDown( PresLeft,  DensLeft,  PresStar );

         GetShockVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensDown_Left, NAN,
			               &ShockVelocity_Left, NULL );

		 VelocityStar = GetVelocityDown( PresLeft, DensLeft, ShockVelocity_Left, PresStar, DensDown_Left );


         RP -> SS.Leftt.Right_Yes       = false;
         RP -> SS.Leftt.ShockVelocity   = ShockVelocity_Left;
         RP -> SS.Leftt.PresUpStream    = PresLeft;
         RP -> SS.Leftt.DensUpStream    = DensLeft;
         RP -> SS.Leftt.VelyUpStream    = VelocityLeft;
         RP -> SS.Leftt.PresDownStream  = PresStar; 
         RP -> SS.Leftt.DensDownStream  = DensDown_Left;
         RP -> SS.Leftt.VelyDownStream  = VelocityStar;


		 DensDown_Right = GetDensDown( PresRight, DensRight, PresStar );

         GetShockVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, NAN, 
		                   NULL, &ShockVelocity_Right );


         RP -> SS.Right.Right_Yes       = true;
         RP -> SS.Right.ShockVelocity   = ShockVelocity_Right;
         RP -> SS.Right.PresUpStream    = PresRight; 
         RP -> SS.Right.DensUpStream    = DensRight;
         RP -> SS.Right.VelyUpStream    = VelocityRight;
         RP -> SS.Right.PresDownStream  = PresStar;
         RP -> SS.Right.DensDownStream  = DensDown_Right;
         RP -> SS.Right.VelyDownStream  = VelocityStar;
	 break;

	 case 2:
         DensDown_Left = GetDensDownRarefaction( PresStar, PresLeft, DensLeft );

		 VelocityStar = GetVelocityDownRarefaction( PresStar, DensDown_Left, PresLeft, DensLeft, VelocityLeft );

	     GetHeadTailVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensDown_Left, VelocityStar,
						      &HeadVelocity_Left, &TailVelocity_Left, false );

         RP -> RS.Leftt.Right_Yes       = false;
         RP -> RS.Leftt.PresUpStream    = PresLeft;
         RP -> RS.Leftt.DensUpStream    = DensLeft;
         RP -> RS.Leftt.VelyUpStream    = VelocityLeft;
         RP -> RS.Leftt.PresDownStream  = PresStar;
         RP -> RS.Leftt.DensDownStream  = DensDown_Left;
         RP -> RS.Leftt.VelyDownStream  = VelocityStar;
         RP -> RS.Leftt.VelocityHead    = HeadVelocity_Left;
         RP -> RS.Leftt.VelocityTail    = TailVelocity_Left;


		 DensDown_Right = GetDensDown( PresRight, DensRight, PresStar );

         GetShockVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, NAN, 
		                   NULL, &ShockVelocity_Right );

		 VelocityStar = GetVelocityDown( PresRight, DensRight, ShockVelocity_Right, PresStar, DensDown_Right );

         RP -> RS.Right.Right_Yes       = true;
         RP -> RS.Right.ShockVelocity   = ShockVelocity_Right;
         RP -> RS.Right.PresUpStream    = PresRight;
         RP -> RS.Right.DensUpStream    = DensRight;
         RP -> RS.Right.VelyUpStream    = VelocityRight;
         RP -> RS.Right.PresDownStream  = PresStar;
         RP -> RS.Right.DensDownStream  = DensDown_Right;
         RP -> RS.Right.VelyDownStream  = VelocityStar;
     break;

	 case 3:
         //DensTail = GetDensInTail( PresStar, PresLeft, DensLeft );
	     //GetHeadTailVelocity( PresLeft, DensLeft, VelocityLeft, PresStar, DensTail, VelocityStar,
		 //   			      &HeadVelocity_Left, &TailVelocity_Left, false );

		 //VelocityStar = GetVelocityInTail( PresStar, DensTail, PresLeft, DensLeft, VelocityLeft );

         //RP -> RS.Leftt.Right_Yes       = false;
         //RP -> RS.Leftt.PresHead        = PresLeft;
         //RP -> RS.Leftt.DensHead        = DensLeft;
         //RP -> RS.Leftt.VelocityHead    = VelocityLeft;
         //RP -> RS.Leftt.PresTail        = PresStar;
         //RP -> RS.Leftt.DensTail        = DensTail;
         //RP -> RS.Leftt.VelocityTail    = VelocityStar;

         //DensTail = GetDensInTail( PresStar, PresRight, DensRight );

         //RP -> RS.Right.Right_Yes       = true;
         //RP -> RS.Right.PresHead        = PresRight;
         //RP -> RS.Right.DensHead        = DensRight;
         //RP -> RS.Right.VelocityHead    = VelocityRight;
         //RP -> RS.Right.PresTail        = PresStar;
         //RP -> RS.Right.DensTail        = DensTail;
         //RP -> RS.Right.VelocityTail    = VelocityStar;
     break;

     default:
		 printf("wave pattern was not found!!\n");
  }


  return Pattern;

}
