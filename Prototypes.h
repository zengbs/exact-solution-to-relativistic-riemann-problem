#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__

int GetWavePattern( struct InitialCondition *IC );

double Velocity_LC ( double PresStar, double DensStarLeft, double PresLeft, double DensLeft, bool Shock );

double Velocity_RC ( double PresStar, double DensStarRight, double PresRight, double DensRight, bool Shock );

double A_Plus ( double Pres, double Dens, double PresLeft, double DensLeft );

double A_Minus ( double Pres, double Dens, double PresRight, double DensRight );

double TaubAdiabatic ( double PresUp, double DensUp, double PresDown );

double PresFunction( double PresStar, struct InitialCondition *IC );

double SoundSpeed( double Pres, double Dens );

double Enthalpy( double Pres, double Dens );

double TotalEngy ( double Pres, double Dens );

double RootFinder( struct InitialCondition *IC, double AbsErr, double RelErr );

#endif
