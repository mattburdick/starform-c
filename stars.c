/*----------------------------------------------------------------------*/
/*                               stars.c                                */
/*                                                                      */
/*  Routines for determining initial star type.  Many of the functions  */
/*  used in this file come from the book "Astrophysics I" by Richard    */
/*  Bowers and Terry Deeming.  Each has been marked with the equation   */
/*  number from that book.                                              */
/*  Another source used is George Abell's "Exploration of the Universe".*/
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
#include	<math.h>
#include    <errno.h>
#include	<stdio.h>

#ifdef MSDOS
#include	<process.h>
#include	<stdlib.h>
#endif

#include	"config.h"
#include	"const.h"
#include	"structs.h"

#ifdef MSDOS
#include	"protos.h"
#endif

extern double about();
extern double power();
extern double random_number();

/*--------------------------------------------------------------------------*/
/*  The 'ms_stardata' array holds information on main-sequence stars (V).   */
/*  Each line contains:                                                     */
/*      {spec_class    spec_num     max_mass	 	percentage}               */
/*                                                                          */
/*  The spectral class and spectral numbers are fairly self-explanatory -   */
/*  the max_mass is the maximum mass a star of that spectral class and      */
/*  number may have.  Mass is a unitless ratio of the star's mass to the    */
/*  Sun's.  Percentage is the percent of all main-sequence stars in the     */
/*  range from that data line to the previous line.  Therefore, 35% of      */
/*  all main-sequence stars in our stellar neighborhood fall in the range   */
/*  from M5 to M9.  Data in these arrays is based on numbers from Bowers'   */
/*  and Deeming's "Astrophysics I", p. 31.                                  */
/*     The data representing percentage distribution of stars for each      */
/*  spectral class is based on a survey of near stars.  Better data should  */
/*  be used some day.                                                       */
/*--------------------------------------------------------------------------*/
spectral_info ms_stardata[] = {
	{'M', 9, 0.1,  0 },
	{'M', 5, 0.2, 35 },
	{'M', 0, 0.5, 36 },
	{'K', 5, 0.7,  7 },
	{'K', 0, 0.8,  7 },
	{'G', 5, 0.9,  3 },
	{'G', 0, 1.1,  3 },
	{'F', 5, 1.3,  2 },
	{'F', 0, 1.7,  1 },
	{'A', 5, 2.0,  1 },
	{'A', 0, 3.2,  1 },
	{'B', 5, 6.5,  1 },
	{'B', 0, 17.8, 1 },
	{'O', 5, 39.8, 1 },
	{'O', 0, 60.0, 1 }
};

/*
 *  The 'wd_stardata' array holds information on white dwarfs (D):
 */
spectral_info wd_stardata[] = {
	{'M', 9, 0.0,  0 },
	{'M', 5, 0.2,  0 },
	{'M', 0, 0.4,  0 },
	{'K', 5, 0.4,  1 },
	{'K', 0, 0.4,  1 },
	{'G', 5, 0.5,  1 },
	{'G', 0, 0.6,  1 },
	{'F', 5, 0.7,  4 },
	{'F', 0, 0.8,  8 },
	{'A', 5, 1.0, 28 },
	{'A', 0, 0.5, 32 },
	{'B', 5, 0.4, 13 },
	{'B', 0, 0.4,  9 },
	{'O', 5, 0.5,  1 },
	{'O', 0, 0.7,  1 }
};

/*
 *  The 'g_stardata' array holds information on giant stars (III):
 */
spectral_info g_stardata[] = {
	{'M', 9, 8.7,  0 },
	{'M', 5, 7.9, 12 },
	{'M', 0, 6.3, 19 },
	{'K', 5, 5.0, 26 },
	{'K', 0, 4.0, 25 },
	{'G', 5, 3.2,  5 },
	{'G', 0, 2.5,  4 },
	{'F', 5, 2.4,  2 },
	{'F', 0, 2.5,  1 },
	{'A', 5, 2.7,  1 },
	{'A', 0, 3.4,  1 },
	{'B', 5, 7.0,  1 },
	{'B', 0, 30.3, 1 },
	{'O', 5, 60.0, 1 },
	{'O', 0, 70.0, 1 }
};

/*
 *  The 'sg_stardata' array holds information on supergiant stars (Ia):
 */
spectral_info sg_stardata[] = {
	{'M', 9, 22.3,  0 },
	{'M', 5, 19.9, 12 },
	{'M', 0, 15.8, 13 },
	{'K', 5, 15.0,  3 },
	{'K', 0, 12.6,  4 },
	{'G', 5, 11.6,  3 },
	{'G', 0, 10.0,  3 },
	{'F', 5, 11.8,  8 },
	{'F', 0, 12.6,  7 },
	{'A', 5, 13.2,  6 },
	{'A', 0, 15.8,  6 },
	{'B', 5, 30.2, 12 },
	{'B', 0, 50.1, 13 },
	{'O', 5, 70.0,  4 },
	{'O', 0, 90.0,  6 }
};

/*--------------------------------------------------------------------------*/
/*   This is eq. 3.52 from "Astrophysics I" by Bowers and Deeming.          */
/*   The mass_ratio is unitless and is a ratio of the stellar mass to that  */
/*   of the Sun.  Both alpha and beta are unitless constants.               */
/*   Note that for a main-sequence G3 star like the Sun, this function      */
/*   overestimates the luminosity slightly.  It does, however, fit the      */
/*   mass-luminosity curve fairly well.                                     */
/*--------------------------------------------------------------------------*/
double luminosity(mass_ratio, lum_class)
double mass_ratio;
int lum_class;
{
	double temp, alpha, beta;

	if ((temp = log10(mass_ratio)) == -HUGE_VAL) {
		perror("luminosity function: stellar mass bad");
		exit(1);
	}
	if (lum_class == MAIN_SEQUENCE) {
/* Set up Bowers and Deeming's alpha and beta constants: */
		if (mass_ratio <=0.5) {
			alpha = 2.85;
			beta = -0.15;
		}
		else if (mass_ratio < 2.5) {
			alpha = 3.6;
			beta = 0.073;
		}
		else {
			alpha = 2.91;
			beta = 0.479;
		}
		temp = beta + alpha * temp;
		temp = power(10.0,temp);
	}
	else if (lum_class == GIANT) {
		temp = temp * 3.3;
		temp = power(10.0, temp);
	}
	else if (lum_class == SUPERGIANT) {
		temp = (temp + 0.22) / 0.33;
		temp = power(10.0, temp);
	}
	else if (lum_class == WHITE_DWARF) {
		temp = mass_ratio * 5.67E-4;
	}
	else temp = 1.0;
	return(temp);
}

/*--------------------------------------------------------------------------*/
/*   This is eq. 3.53 from "Astrophysics I" by Bowers and Deeming.          */
/*   The mass_ratio is unitless and is a ratio of the stellar mass to that  */
/*   of the Sun.  The stellar radius returned is in units of AU.            */
/*--------------------------------------------------------------------------*/
double star_radius (mass_ratio, lum_class, cooler_than_G0)
double mass_ratio;
int lum_class, cooler_than_G0;
{
	double temp;

	if ((temp = log10(mass_ratio)) == -HUGE_VAL) {
		perror("stellar radius function: stellar mass bad");
		exit(1);
	}
	if (lum_class == MAIN_SEQUENCE) {
		if (mass_ratio <= 0.4) {
			temp = temp + 0.1;
			temp = power(10.0, temp);
		}
		else {
			temp = 0.73 * temp;
			temp = power(10.0, temp);
		}
	}
	else if (lum_class == GIANT) {
		temp = temp * 2.0;
		temp = power(10.0, temp);
	}
	else if (lum_class == SUPERGIANT) {
		if (cooler_than_G0) {
			temp = (temp - 0.32) / 0.34;
			temp = power(10.0, temp);
		}
		else {
			temp = (temp - 2.7) / -0.86;
			temp = power(10.0, temp);
		}
	}
	else if (lum_class == WHITE_DWARF) {
		temp = about(0.02, 0.005);
	}
	else temp = 1.0;
	temp = temp * SOLAR_RADII_PER_AU;  /* convert to AU units */
	return(temp);
}

/*--------------------------------------------------------------------------*/
/*   Both the main sequence lifetime and the age returned are in units of   */
/*   years.  The lifetime passed in is guaranteed to be >= 1 million.       */
/*--------------------------------------------------------------------------*/
double star_age (lifetime)
double lifetime;
{
	double temp;

	if (lifetime >= 6.0E9)
		temp = random_number(1.0E9, 6.0E9);
	else if (lifetime > 1.0E9)
		temp = random_number(1.0E9, lifetime);
	else temp = random_number(1.0E6, lifetime);
	return(temp);
}

/*--------------------------------------------------------------------------*/
/*   Using the information in the 'stardata' array, we can determine        */
/*   what spectral class and spectral number to apply to a star (given the  */
/*   stellar mass ratio of that star).  This function searches 'stardata'   */
/*   for the correct spectral class catagory, then calculates the spectral  */
/*   number.                                                                */
/*--------------------------------------------------------------------------*/
char *classify (mass_ratio, lum_class)
double mass_ratio;
int lum_class;
{
	spectral_info *stardata;
	int i, modifier, temp;
	double prev_mass;
	char buf[CLASSIFICATION_SIZE];

	switch (lum_class) {
		case GIANT:
			stardata = g_stardata;
			break;
		case SUPERGIANT:
			stardata = sg_stardata;
			break;
		case WHITE_DWARF:
			stardata = wd_stardata;
			break;
		case MAIN_SEQUENCE:
		default:
			stardata = ms_stardata;
			break;
	}
	for (i = 0, prev_mass = 0.049;
		((i < 14) && (stardata[i].max_mass < mass_ratio));
		prev_mass = stardata[i].max_mass, i++)
	;
	if ((i == 14) && (stardata[i].max_mass < mass_ratio)) {
		sprintf(buf, "%c%c %c\0", '?', '?', '?');
	}
	else {
		temp = (int) (5.0 * (stardata[i].max_mass - mass_ratio) /
				   (stardata[i].max_mass - prev_mass));
		modifier = stardata[i].spec_num + temp;
		switch (lum_class) {
			case GIANT:
				sprintf(buf, "%c%d III\0", stardata[i].spec_class, modifier);
				break;
			case SUPERGIANT:
				sprintf(buf, "%c%d Ia\0", stardata[i].spec_class, modifier);
				break;
			case WHITE_DWARF:
				sprintf(buf, "%c%d D\0", stardata[i].spec_class, modifier);
				break;
			case MAIN_SEQUENCE:
			default:
				sprintf(buf, "%c%d V\0", stardata[i].spec_class, modifier);
				break;
		}
	}
	return(buf);
}

/*--------------------------------------------------------------------------*/
/*   Given the luminosity class (giant, white dwarf, etc), spectral class,  */
/*   and spectral number of a star, this function will return the mass      */
/*   interpolated from the stardata structure arrays.                       */
/*--------------------------------------------------------------------------*/
double star_mass(lum_class, spec_class, spec_num)
char spec_class;
int spec_num, lum_class;
{
	spectral_info *stardata;
	int i;
	double prev_mass, temp;

	switch (lum_class) {
		case GIANT:
			stardata = g_stardata;
			break;
		case SUPERGIANT:
			stardata = sg_stardata;
			break;
		case WHITE_DWARF:
			stardata = wd_stardata;
			break;
		case MAIN_SEQUENCE:
		default:
			stardata = ms_stardata;
			break;
	}
	for (i = 0, prev_mass = 0.0;
		((i < 14) && (stardata[i].spec_class != spec_class)
		|| (stardata[i].spec_num > spec_num));
		prev_mass = stardata[i].max_mass, i++)
	;
	if ((i == 14) && (stardata[i].spec_class != spec_class)) {
		return(0.0);  /* Must have been an error somewhere */
	}
	else {
		temp = ((spec_num - stardata[i].spec_num)
				* (stardata[i].max_mass - prev_mass)) / 5.0;
		temp = stardata[i].max_mass - temp;
		return(temp);
	}
}

/*--------------------------------------------------------------------------*/
/*   This function checks the information provided by the user via the      */
/*   '-t' command-line flag.  An error indicator will be returned if the    */
/*   command-line info is incorrect.                                        */
/*--------------------------------------------------------------------------*/
int verify_startype (lum_id, spec_num, spec_class)
char lum_id, spec_class;
int spec_num;
{
	int error_type = 0;

	if ((spec_class != 'O') && (spec_class != 'B')
	    && (spec_class != 'A') && (spec_class != 'F')
	    && (spec_class != 'G') && (spec_class != 'K')
	    && (spec_class != 'M')) {
	    error_type |= BAD_SPECTRA;
	}
	if ((spec_num < 0) || (spec_num > 9)) {
		error_type |= BAD_MOD;
	}
	if ((lum_id != 'M') && (lum_id != 'G')
	    && (lum_id != 'S') && (lum_id != 'D')) {
	    error_type |= BAD_LUMINOSITY;
	}
	return(error_type);
}

/*--------------------------------------------------------------------------*/
/*   The startype_error function interprets and displays an error           */
/*   indicator returned by the verify_startype function.                    */
/*--------------------------------------------------------------------------*/
void startype_error(errornum, spec_class, spec_num, lum_id)
int errornum, spec_num;
char spec_class, lum_id;
{
	if (errornum & BAD_SPECTRA)
		printf("ERROR: invalid spectral class <%c>\n", spec_class);
	if (errornum & BAD_MOD)
		printf("ERROR: invalid spectral class modifier <%d>\n", spec_num);
	if (errornum & BAD_LUMINOSITY)
		printf("ERROR: invalid luminosity class <%c>\n", lum_id);
}

/*--------------------------------------------------------------------------*/
/*   This function uses the stardata array passed to it to determine the    */
/*   mass of a random star.                                                 */
/*   Mass is returned as a ratio of the star's mass to the Sun's.           */
/*--------------------------------------------------------------------------*/
double rand_star_mass (startype)
int startype;
{
	spectral_info *stardata;
	int temp, percent = 0, i;
	double mass, prev_mass = 0.0;

	switch (startype) {
		case GIANT:
			stardata = g_stardata;
			break;
		case SUPERGIANT:
			stardata = sg_stardata;
			break;
		case WHITE_DWARF:
			stardata = wd_stardata;
			break;
		case MAIN_SEQUENCE:
		default:
			stardata = ms_stardata;
			break;
	}
	temp = (int)random_number(0.0, 100.0);
	for (i=0; (i <= 14); i++) {
		percent = percent + stardata[i].percentage;
		if (temp <= percent) {
			mass = random_number(stardata[i].max_mass, prev_mass);
			return(mass);
		}
		prev_mass = stardata[i].max_mass;
	}
	return (1.0);  /* should never get here */
}

/*--------------------------------------------------------------------------*/
/*   According to George Abell's "Exploration of the Universe", (fourth     */
/*   edition), about 90% of all stars in the local neighborhood are main-   */
/*   sequence stars, while about 10% are white dwarfs and less than 1% are  */
/*   giants or supergiants.  This function reflects those percentages.  If  */
/*   you are interested in larger stars, you can always generate them       */
/*   using the '-t' flag!                                                   */
/*--------------------------------------------------------------------------*/
int rand_type () {
	int temp;

	temp = (int)random_number(0.0, 100.0);
	if (temp <= 1) {                           /* Giant or supergiant  */
		temp = (int)random_number(0.0, 100.0);
		if (temp <=70) {                      /* Must be a giant      */
			return(GIANT);
		}
		else {						   /* Must be a supergiant */
			return(SUPERGIANT);
		}
	}
	else if (temp <= 10) {                     /* White dwarf          */
		return(WHITE_DWARF);
	}
	else {                                     /* Main sequence        */
		return(MAIN_SEQUENCE);
	}
}
