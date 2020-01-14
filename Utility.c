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

  PresStar = RootFinder( PresFunction, (void*)IC, 0.0, __DBL_EPSILON__, 0.5, 1e-2, 1.0 );


  double ShockVelocity_Left,  DensDown_Left;
  double ShockVelocity_Right, DensDown_Right;

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



//         RP -> RS.Right.Right_Yes       = true;
//         RP -> RS.Right.PresHead        = 
//         RP -> RS.Right.DensHead        = 
//         RP -> RS.Right.VelocityHead    = 
//         RP -> RS.Right.PresTail        = 
//         RP -> RS.Right.DensTail        = 
//         RP -> RS.Right.VelocityTail    = 
//         RP -> RS.Right.Xi              = 
//
//
//		 DensDown_Right = GetDensDown( PresRight, DensRight, PresStar );
//
//         GetShockVelocity( PresRight, DensRight, VelocityRight, PresStar, DensDown_Right, NAN, 
//						NULL, &ShockVelocity_Right );
//
//
//         RP -> RS.Leftt.Right_Yes       = false;
//         RP -> RS.Leftt.ShockVelocity   = 
//         RP -> RS.Leftt.PresUpStream    = 
//         RP -> RS.Leftt.DensUpStream    = 
//         RP -> RS.Leftt.VelyUpStream    = 
//         RP -> RS.Leftt.PresDownStream  = 
//         RP -> RS.Leftt.DensDownStream  = 
//         RP -> RS.Leftt.VelyDownStream  = 
    	 break;

	 case 3:
//         RP -> RR.Right.Right_Yes       = 
//         RP -> RR.Right.PresHead        = 
//         RP -> RR.Right.DensHead        = 
//         RP -> RR.Right.VelocityHead    = 
//         RP -> RR.Right.PresTail        = 
//         RP -> RR.Right.DensTail        = 
//         RP -> RR.Right.VelocityTail    = 
//         RP -> RR.Right.Xi              = 
//
//         RP -> RR.Leftt.Right_Yes       = 
//         RP -> RR.Leftt.PresHead        = 
//         RP -> RR.Leftt.DensHead        = 
//         RP -> RR.Leftt.VelocityHead    = 
//         RP -> RR.Leftt.PresTail        = 
//         RP -> RR.Leftt.DensTail        = 
//         RP -> RR.Leftt.VelocityTail    = 
//         RP -> RR.Leftt.Xi              = 
	     break;
     default:
		 printf("wave pattern was not found!!\n");
  }


  return Pattern;

}
