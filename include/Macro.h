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


// shock regions
#define N_REGIONS  6
#define L_ORIGINAL 1    // original left state
#define L_FAN      2    // rarefaction fan on the left
#define L_CONTACT  3    // left region of contact
#define R_CONTACT  4    // right region of contact
#define R_FAN      5    // rarefaction fan on the right
#define R_ORIGINAL 6    // original right state

#endif // #ifndef __MACRO_H__
