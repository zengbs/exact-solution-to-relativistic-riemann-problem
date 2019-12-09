#ifndef __STRUCT_H__
#define __STRUCT_H__

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

struct RareFaction
{
  bool   Right_Yes;
  double PresHead;
  double DensHead;
  double VelocityHead;
  double PresTail;
  double DensTail;
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
  struct RareFaction Right;
  struct Shock       Leftt;
};

struct RRWaves
{
  struct RareFaction Right;
  struct RareFaction Leftt;
};

struct RiemannProblem
{
  struct SSWaves SS;
  struct RSWaves RS;
  struct RRWaves RR;
};

#endif
