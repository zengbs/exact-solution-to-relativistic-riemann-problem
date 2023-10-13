#ifndef __MACRO_H__
#define __MACRO_H__

// funcitons
#define SQR( x )        ( ( x ) * ( x ) )
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )
#define SIGN( a )       (  ( (a) < 0.0 ) ? -1.0 : +1.0  )

// EoS types
#define TM    0
#define GAMMA 1

#endif // #ifndef __MACRO_H__
