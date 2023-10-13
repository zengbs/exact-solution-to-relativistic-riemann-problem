#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "Struct.h"
#include "Prototypes.h"

double RootFinder( double(*Function)(double X, void *params) , void *params, double AbsErr, double RelErr,
			       double Guess, double LowerBound, double UpperBound, const char FunctionName[] )
{
// Make sure the root is between LowerBound and UpperBound
  while ( Function(LowerBound, params) *  Function(UpperBound, params) >= 0.0 )
  {
    if ( fabs(Function(LowerBound, params)) > fabs(Function(UpperBound, params)) )
    {
       UpperBound *= 10.0;
    }
    else
    {
       LowerBound /= 10.0;
    }

    if ( Function(LowerBound, params) == Function(UpperBound, params) )
    {
       UpperBound *= 10.0;
       LowerBound *= 10.0;
    }
    //printf("LowerBound=%e, UpperBound=%e, f(LowerBound)=%e, f=(UpperBound)=%e: %s\n",
    //LowerBound, UpperBound, Function(LowerBound, params), Function(UpperBound, params), FunctionName );
  }
//printf("hi\n");
  int status;

  int iter = 0, max_iter = 100;

  const gsl_root_fsolver_type *T;

  gsl_root_fsolver *s;

  double Root, RootTemp;

  gsl_function F;

  F.function = Function;

  F.params = params;

  T = gsl_root_fsolver_brent;

  s = gsl_root_fsolver_alloc (T);

  gsl_root_fsolver_set (s, &F, LowerBound, UpperBound );

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

  gsl_root_fsolver_free(s);

  return Root;
}
