/*----------------------------------------------------------------------*/
/*                               utils.c                                */
/*                                                                      */
/*  General-purpose utility routines.  Note that if your machine is     */
/*  having difficulty is the random-number generation routine in some   */
/*  way, the problem will often manifest itself here in the form of     */
/*  the "error in power" message.  If so, try using a debugger to       */
/*  trace the 'random_number' function to see what is passed in and     */
/*  returned.                                                           */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
#include    <math.h>
#include    <errno.h>
#include	<stdio.h>
#include    <stdlib.h>

#ifdef MSDOS
#include	<stdlib.h>
#endif

#include    "config.h"
#include    "const.h"
#include    "structs.h"

#ifdef MSDOS
#include    "protos.h"
#endif

/* extern int errno;*/

double power(x, y)
double x,y;
{
    double result;

    errno = 0;
    if (((result = pow(x,y)) == 0.0) && (errno == EDOM)) {
        printf("error in power\n");
    }
    return (result);
}

/*----------------------------------------------------------------------*/
/*  This function returns a random real number between the specified    */
/*  real number bounds.                                                 */
/*----------------------------------------------------------------------*/

double random_number(bound1, bound2)
double bound1, bound2;
{
    double range, lowbound;

    if (bound1 > bound2) {
        range = bound1 - bound2;
        lowbound = bound2;
    }
    else if (bound2 > bound1) {
        range = bound2 - bound1;
        lowbound = bound1;
    }
    else {
        return(bound1);    /* Since bound1 must equal bound2 */
    }
    return((((double)rand()) / (double)(RAND_MAX)) * range + lowbound);
}


/*----------------------------------------------------------------------*/
/*   This function returns a value within a certain variation of the    */
/*   exact value given it in 'value'.                                   */
/*----------------------------------------------------------------------*/

double about(value, variation)
double value, variation;
{
    return(value + (value * random_number(-variation,variation)));
}

double random_eccentricity()
{
    return(1.0 - power(random_number(0.0001, 1.0),ECCENTRICITY_COEFF));
}


