#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "Struct.h"
#include "Prototypes.h"

double RootFinder( struct InitialCondition *IC, double AbsErr, double RelErr )
{
  double DensLeft      = IC -> DensLeft     ;
  double VelocityLeft  = IC -> VelocityLeft ;
  double PresLeft      = IC -> PresLeft     ;
  double DensRight     = IC -> DensRight    ;
  double VelocityRight = IC -> VelocityRight;
  double PresRight     = IC -> PresRight    ;

  int status;

  int iter = 0, max_iter = 50;

  const gsl_root_fsolver_type *T;

  gsl_root_fsolver *s;

  double Guess = 0.5 * ( PresLeft + PresRight );

  double Root, RootTemp;

  gsl_function F;

  F.function = &PresFunction;

  F.params = IC;

  T = gsl_root_fsolver_brent;

  s = gsl_root_fsolver_alloc (T); 

  gsl_root_fsolver_set (s, &F, PresRight+0.0000000001, PresLeft-0.00000001);

  Root = Guess;

  do
  {
    iter++;

    status = gsl_root_fsolver_iterate (s);

	RootTemp = Root;

    Root = gsl_root_fsolver_root (s);

    status = gsl_root_test_delta (Root, RootTemp, AbsErr, RelErr );
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free;
  
  return Root;
}
