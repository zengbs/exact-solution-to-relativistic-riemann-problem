#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__

#include "struct.h"


void Load_Parameter();

int GetWavePattern( struct InitialCondition *IC );

double JumpConditionForEnthalpy( double EnthalpyDown, void* params );

double Velocity_LC ( double PresStar, double DensStarLeft, double PresLeft, double DensLeft, double VelocityLeft, bool Shock );

double Velocity_RC ( double PresStar, double DensStarRight, double PresRight, double DensRight, double VelocityRight, bool Shock );

double GetEnthalpyDown ( double PresUp, double DensUp, double PresDown );

double PresFunction( double PresStar, void * );

double Flu_SoundSpeed( double Temp );

double Flu_Enthalpy( double Pres, double Dens );

double Flu_TotalInternalEngy ( double Pres, double Dens );

double Enthalpy2Temperature( double Enthalpy );

double RootFinder( double(*Function)(double X, void *params) , void *params, double AbsErr, double RelErr,
                   double Guess, double LowerBound, double UpperBound, const char s[] );


double MassCurrent( double PresUp, double DensUp, double PresDown, double DensDown );

int GetAllInfomation( struct InitialCondition *, struct RiemannProblem * );

void Plot( int Pattern, struct RiemannProblem *, struct PlotParams );

double GetVelocityDown( double PresUp,   double DensUp, double ShockFrontVelocity,
                        double PresDown, double DensDown );


double GetDensDown( double PresUp, double DensUp, double PresDown );

void GetShockVelocity( double PresUp, double DensUp, double V_Up, double PresDown, double DensDown, double *Vs_Left, double *Vs_Right );

void QuadraticSolver( double A, double B, double C , double *PlusRoot, double *MinusRoot );

void GetHeadTailVelocity( double PresHead, double DensHead, double VelocityHead, double PresTail, double DensTail, double VelocityTail,
                          double *HeadVelocity, double *TailVelocity, bool Right_Yes );

double GetDensDownRarefaction( double PresDown, double PresUp, double DensUp );

double GetVelocityDownRarefaction( double PresDown, double DensDown, double PresUp, double DensUp, double VelocityUp, bool Right_Yes );

double FanFunction ( double Dens_at_Xi, void *params );

double GetDensInFan( struct Rarefaction *Rarefaction );

double GetPresInFan( double Dens_at_Xi, double PresUp, double DensUp );

double GetVelocityInFan( double Xi,  double Dens_at_Xi, double Pres_at_Xi, bool Right_Yes );

double U2V( double U );

double Isentropic_Constant ( double Init_Temp, double Init_Dens );

double Isentropic_Dens2Temperature ( double Dens, double Init_Temp, double Init_Dens );

double Isentropic_Temperature2Dens ( double Temperature, double Init_Temp, double Init_Dens );

double Isentropic_Pres2Temperature ( struct Rarefaction *Rarefaction );

double Isentropic_Temperature2Pres ( double Temperature, void *params );

double Isentropic_Pres2Dens ( struct Rarefaction *Rarefaction );

double Isentropic_Dens2Pres ( double Dens, double Init_Temp, double Init_Dens );

void RelativeVelocity( double Ua, double Ub, double *LorentzFactor_ab, double *U_ab );

double Isentropic_Dens2Velocity ( double DensDown, struct Rarefaction *upstream );

#endif // #ifndef __PROTOTYPE_H__
