#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__

#include "Struct.h"

int GetWavePattern( struct InitialCondition *IC );

double Velocity_LC ( double PresStar, double DensStarLeft, double PresLeft, double DensLeft, bool Shock );

double Velocity_RC ( double PresStar, double DensStarRight, double PresRight, double DensRight, bool Shock );

double A_PlusFun ( double Temp );

double A_MinusFun ( double Temp);

double TaubAdiabatic ( double PresUp, double DensUp, double PresDown );

double PresFunction( double PresStar, void * );

double Flu_SoundSpeed( double Temp );

double Flu_Enthalpy( double Pres, double Dens );

double Flu_TotalEngy ( double Pres, double Dens );

double RootFinder( double(*Function)(double X, void *params) , void *params, double AbsErr, double RelErr,
			       double Guess, double LowerBound, double UpperBound );


double MassCurrent( double PresUp, double DensUp, double PresDown, double DensDown );

int GetAllInfomation( struct InitialCondition *, struct RiemannProblem  * );
void Plot( int Pattern, struct RiemannProblem *, struct PlotParams   );

double GetVelocityDown( double PresUp,   double DensUp, double ShockFrontVelocity,
                        double PresDown, double DensDown );


double GetDensDown( double PresUp, double DensUp, double PresDown  );

void GetShockVelocity( double PresUp,   double DensUp,   double V_Up, 
				    double PresDown, double DensDown,
			        double *Vs_Left, double *Vs_Right );

void QuadraticSolver( double A, double B, double C , double *PlusRoot, double *MinusRoot);

void GetHeadTailVelocity( double PresHead, double DensHead, double VelocityHead,
			              double PresTail, double DensTail, double VelocityTail,
                          double *HeadVelocity, double *TailVelocity, bool Right_Yes );

double GetDensDownRarefaction( double PresDown, double PresUp, double DensUp );

double GetVelocityDownRarefaction( double PresDown, double DensDown, double PresUp, double DensUp, double VelocityUp, bool Right_Yes );

double GetDensInFan( double Cs2, double PresHead, double DensHead );


double GetPresInFan( double DensInFan, double PresHead, double DensHead );

double GetVelocityInFan( double Cs2, double Xi, bool Right_Yes );

double GetSoundSpeedInFan ( struct Rarefaction *Rarefaction );

double TemperatureFunction ( double Temp, void *params );

double U2V(double U);

#endif
