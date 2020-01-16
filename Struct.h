#ifndef __STRUCT_H__
#define __STRUCT_H__

#include <stdbool.h>

struct InitialCondition
{
  double DensLeft;
  double VelocityLeft;
  double PresLeft;
  double DensRight;
  double VelocityRight;
  double PresRight;
};


struct Shock
{
  bool   Right_Yes;
  double ShockVelocity;
  double PresUpStream;
  double DensUpStream;
  double VelyUpStream;
  double PresDownStream;
  double DensDownStream;
  double VelyDownStream;
};

struct Rarefaction
{
  bool   Right_Yes;
  double PresUpStream;
  double DensUpStream;
  double VelyUpStream;
  double PresDownStream;
  double DensDownStream;
  double VelyDownStream;
  double VelocityHead;
  double VelocityTail;
  double Xi;
};

struct SSWaves
{
  struct Shock Right;
  struct Shock Leftt;
};

struct RSWaves
{
  struct Rarefaction Leftt;
  struct Shock       Right;
};

struct RRWaves
{
  struct Rarefaction Right;
  struct Rarefaction Leftt;
};

struct RiemannProblem
{
  struct SSWaves SS;
  struct RSWaves RS;
  struct RRWaves RR;
};

struct PlotParams
{
  double DT;
  double End_T;
  double X_Left;
  double X_Right;
  int NCell;
};

#endif
