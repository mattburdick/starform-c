/*----------------------------------------------------------------------*/
/*                               config.h                               */
/*                                                                      */
/*  Header file for resolving machine-dependent parameters such as the  */
/*  type of 'rand' function your machine provides.                      */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/

#define VERSION "3.6"


/*
 * For machines on which the rand() function returns a long int (32 bits).
 * Suns, for instance:
 */

#if defined(sun) || defined(NeXT)
#define LONG_RAND
#endif


/*
 * For Berkeley-based C compilers.  Many Berkeley-ish machines, such as Suns,
 * HP-UX machines, and Apollos are compatible enough with SYSV that this flag
 * is not needed.  If your machine has '/usr/include/strings.h' rather than
 * 'string.h', you need this defined.
 */

/*#define BSD*/

