#ifndef __MACRO_H__
#define __MACRO_H__

#include <stdio.h>

#define LOCATION __FILE__, __FUNCTION__, __LINE__
#define REPORT_ERROR { printf( "Error: %s(%s):%d\n", LOCATION ); abort(); }


#define SQR( x ) ( ( x ) * ( x ) )
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )
#define SIGN( a )       (  ( (a) < 0.0 ) ? -1.0 : +1.0  )
#define TM    0
#define GAMMA 1

#endif
