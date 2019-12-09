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


struct RareFactionFan
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

#endif
