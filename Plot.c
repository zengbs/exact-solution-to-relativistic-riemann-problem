#include <stdio.h>
#include "Struct.h"
#include "Prototypes.h"


void Plot( int Pattern, struct RiemannProblem *RP, struct PlotParams plot )
{

  FILE *fptr[999999];

  double DT      = plot.DT     ;
  double End_T   = plot.End_T  ;
  double X_Left  = plot.X_Left ;
  double X_Right = plot.X_Right;
  int    NCell   = plot.NCell  ;

  double dX = (X_Left-X_Right)/(double)NCell;
  double X0 = 0.5*( X_Left + X_Right );

  char fileName[100];
  sprintf (fileName, "./000000.dat");
  fptr[0] = fopen (fileName, "w");

  double time = 0.0;
  fprintf (fptr[0], "# time = %e\n", time);
  fprintf (fptr[0], "# x[1]       d1[2]       v1[3]       p1[4]\n");

  fprintf (fptr[0], "%e %e %e %e\n", X_Left,  RP->SS.Leftt.DensUpStream, RP->SS.Leftt.VelyUpStream, RP->SS.Leftt.PresUpStream);
  fprintf (fptr[0], "%e %e %e %e\n", X0-dX,   RP->SS.Leftt.DensUpStream, RP->SS.Leftt.VelyUpStream, RP->SS.Leftt.PresUpStream);

  fprintf (fptr[0], "%e %e %e %e\n", X0,      RP->SS.Leftt.DensUpStream, RP->SS.Leftt.VelyUpStream, RP->SS.Leftt.PresUpStream);
  fprintf (fptr[0], "%e %e %e %e\n", X_Right, RP->SS.Leftt.DensUpStream, RP->SS.Leftt.VelyUpStream, RP->SS.Leftt.PresUpStream);



  int j = 1;

  double X_Max, X_Min;
  int Ng_1, Ng_2, Ng_3, Ng_4;

  time = time + DT;

  while (time <= End_T)
  {
    sprintf (fileName, "./%06d.dat", j);
    fptr[j] = fopen (fileName, "w");
    fprintf (fptr[j], "# time = %e\n", time);
    fprintf (fptr[j], "# x=[1], d1=[2], v1=[3], p1=[4], D1=[5], M1=[6], E1=[7], Lorentz fac.=[8]\n");

    int region;

    for ( region = 1; region <= 4; region++)
    {
        switch (region)
        {
          case 1:
          X_Max = RP->SS.Leftt.ShockVelocity * time + X0;
          X_Min = 0.0;

          fprintf (fptr[j], "%e %e %e %e\n", X_Min,    RP->SS.Leftt.DensUpStream, RP->SS.Leftt.VelyUpStream, RP->SS.Leftt.PresUpStream );
          fprintf (fptr[j], "%e %e %e %e\n", X_Max-dX, RP->SS.Leftt.DensUpStream, RP->SS.Leftt.VelyUpStream, RP->SS.Leftt.PresUpStream );
          break;

          case 2:
          X_Max = RP->SS.Leftt.VelyDownStream * time + X0;
          X_Min = RP->SS.Leftt.ShockVelocity  * time + X0;

          fprintf (fptr[j], "%e %e %e %e\n", X_Min,   RP->SS.Leftt.DensDownStream, RP->SS.Leftt.VelyDownStream, RP->SS.Leftt.PresDownStream );
          fprintf (fptr[j], "%e %e %e %e\n", X_Max-dX,RP->SS.Leftt.DensDownStream, RP->SS.Leftt.VelyDownStream, RP->SS.Leftt.PresDownStream );
          break;

          case 3:
          X_Max = RP->SS.Right.ShockVelocity * time + X0;
          X_Min = RP->SS.Right.VelyDownStream* time + X0;

          fprintf (fptr[j], "%e %e %e %e\n", X_Min,   RP->SS.Right.DensDownStream, RP->SS.Right.VelyDownStream, RP->SS.Right.PresDownStream );
          fprintf (fptr[j], "%e %e %e %e\n", X_Max-dX,RP->SS.Right.DensDownStream, RP->SS.Right.VelyDownStream, RP->SS.Right.PresDownStream );
          break;

          case 4:
          X_Max =  X_Right;
          X_Min =  RP->SS.Right.ShockVelocity * time + X0;

          fprintf (fptr[j], "%e %e %e %e\n", X_Min,    RP->SS.Right.DensUpStream, RP->SS.Right.VelyUpStream, RP->SS.Right.PresUpStream );
          fprintf (fptr[j], "%e %e %e %e\n", X_Max-dX, RP->SS.Right.DensUpStream, RP->SS.Right.VelyUpStream, RP->SS.Right.PresUpStream);
          break;

        }
    }

    fclose (fptr[j]);

    j++;

    time = j * DT;
  }




  printf("Leftt.ShockVelocity  %e\n", RP->SS.Leftt.ShockVelocity  );
  printf("Leftt.PresUpStream   %e\n", RP->SS.Leftt.PresUpStream   ); 
  printf("Leftt.DensUpStream   %e\n", RP->SS.Leftt.DensUpStream   ); 
  printf("Leftt.VelyUpStream   %e\n", RP->SS.Leftt.VelyUpStream   ); 
  printf("Leftt.PresDownStream %e\n", RP->SS.Leftt.PresDownStream ); 
  printf("Leftt.DensDownStream %e\n", RP->SS.Leftt.DensDownStream ); 
  printf("Leftt.VelyDownStream %e\n", RP->SS.Leftt.VelyDownStream ); 
                               
                               
  printf("Right.ShockVelocity  %e\n", RP->SS.Right.ShockVelocity  );
  printf("Right.PresUpStream   %e\n", RP->SS.Right.PresUpStream   ); 
  printf("Right.DensUpStream   %e\n", RP->SS.Right.DensUpStream   ); 
  printf("Right.VelyUpStream   %e\n", RP->SS.Right.VelyUpStream   ); 
  printf("Right.PresDownStream %e\n", RP->SS.Right.PresDownStream ); 
  printf("Right.DensDownStream %e\n", RP->SS.Right.DensDownStream ); 
  printf("Right.VelyDownStream %e\n", RP->SS.Right.VelyDownStream ); 


}
