/*----------------------------------------------------------------------*/
/*                              display.c                               */
/*                                                                      */
/*  Output routines - replace this file if you want to write a          */
/*  graphical interface.                                                */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
#include	<stdio.h>
#include    <math.h>

#ifdef BSD
#include	<strings.h>
#else
#include	<string.h>
#endif

#ifdef MSDOS
#include    <stddef.h>
#include    <malloc.h>
#include	<stdlib.h>
#include    <float.h>
#endif

#include	"config.h"
#include	"const.h"
#include	"structs.h"

#ifdef MSDOS
#include	"protos.h"
#endif

extern double power();
extern double random_number();
extern int flag_graphic;
extern int flag_moons;
extern int flag_tec;

void draw_system(sys)
sys_pointer sys;
{
}

void describe_star (star)
star_pointer star;
{
	printf("Stellar Classification:      %7s\n", star->star_type);
	printf("Stellar mass:                %7.2f solar masses\n",
	       star->stell_mass_ratio);
	printf("Stellar radius:              %7.4f AU\n", star->stell_radius);
	printf("Stellar luminosity:          %7.3f\n", star->stell_luminosity_ratio);
	printf("Age:                         %7.3f billion years\n", (star->age /1.0E9));
	if (star->lum_type == MAIN_SEQUENCE)
		printf("Years left on Main Sequence: %7.3f billion years\n",
			  (star->main_seq_life - star->age) / 1.0E9);
	printf("Earthlike insolation at:     %7.3f AU\n",star->r_ecosphere);
}

void describe_system(sys)
sys_pointer sys;
{
	planet_pointer node1;
	planet_pointer node2;
	star_pointer star_ptr;
	int counter1, counter2;

	printf("                         SYSTEM  CHARACTERISTICS\n\n");
	printf("        PRIMARY STAR\n");
	describe_star(sys->primary_star);
	if (sys->primary_star->next_star != NULL) {
		printf("\n");
		printf("Companion stars present at:\n");
		for (star_ptr=sys->primary_star->next_star, counter1=1;
				star_ptr != NULL;
				star_ptr=star_ptr->next_star) {
			printf("%d\t%7.3lf \t AU\n", counter1, star_ptr->orbit_radius);
			counter1++;
		}
	}
	printf("\n");
	printf("Planets present at:\n");
	for (node1=sys->inner_planet, counter1=1;
	     node1 != NULL;
	     node1=node1->next_planet) {
		if (node1->mass_type != STAR) {
            if (node1->mass_type == GAS_GIANT) {
			    printf("%d\t%7.3lf \t AU  * Gas giant *\n", counter1,
                                                            node1->a);
            }
            else {
			    printf("%d\t%7.3lf \t AU\n", counter1, node1->a);
            }
			counter1++;
		}
	}
	printf("\n\n\n");
	/*
	 *  Loop through the planets, displaying each.  Start with the second
	 *  planet since the first 'planet' is really the primary star.
	 */
	for (node1=sys->inner_planet->next_planet, counter1=1;
	     node1 != NULL;
	     node1=node1->next_planet) {
		if (node1->mass_type == STAR) {
			printf("COMPANION STAR\n");
			printf("Orbital Radius:           %9.3f AU\n",
		       node1->star_ptr->orbit_radius);
			describe_star(node1->star_ptr);
			printf("\n\n");
			continue;						/* skip to next planet */
		}
		printf("Planet %d",counter1++);
		/*
		 *	Continue with this only if we're talking about a gas giant
		 *	or normal planet.
		 */
		if (node1->mass_type == GAS_GIANT) {
			printf("\t*gas giant*\n");
		}
		else printf("\n");
	    if ((int)node1->day == (int)(node1->orb_period * 24.0))
	     	printf("Planet tidally locked (one face to star).\n");
	    if (node1->resonant_period)
	     	printf("Planet almost tidally locked with star\n");
		printf("   Orbital Radius:           %9.3f AU\n",
		       node1->a);
		printf("   Mass:                     %9.3f Earth masses\n",
		       node1->mass * SUN_MASS_IN_EARTH_MASSES);
		if (node1->mass_type == PLANET) {
			printf("   Surface gravity:          %9.2f Earth gees\n",
			       node1->surf_grav);
			printf("   Surface pressure:         %9.3f Earth atm",
			       (node1->surf_pressure / 1000.0));
			if ((node1->greenhouse_effect)
			    && (node1->surf_pressure > 0.0))
				printf("   GREENHOUSE EFFECT\n");
			else printf("\n");
			printf("   Surface temperature:      %9.2f deg Cel\n",
			       (node1->surf_temp -KELVIN_CELCIUS_DIFFERENCE));
		}
	    printf("   Equatorial radius:        %9.1f Km\n",node1->radius);
		printf("   Density:                  %9.3f grams/cc\n",node1->density);
		printf("   Eccentricity of orbit:    %9.3f\n",node1->e);
		printf("   Escape Velocity:          %9.2f Km/sec\n",
		    node1->esc_velocity / CM_PER_KM);
		printf("   Molecular weight retained:%9.2f and above\n",
		    node1->molec_weight);
		printf("   Surface acceleration:     %9.2f cm/sec2\n",
		    node1->surf_accel);
		printf("   Axial tilt:               %9d degrees\n",node1->axial_tilt);
		printf("   Planetary albedo:         %9.3f\n",node1->albedo);
		printf("   Length of year:           %9.2f days\n",
		    node1->orb_period);
		printf("   Length of day:            %9.2f hours\n",node1->day);
		if (node1->mass_type == PLANET) {
			printf("   Boiling pt. of water:     %9.1f deg Cel\n",
			       (node1->boil_point-KELVIN_CELCIUS_DIFFERENCE));
			printf("   Hydrosphere percentage:   %9.2f\n",
			       (node1->hydrosphere * 100.0));
			printf("   Cloud cover percentage:   %9.2f\n",
			       (node1->cloud_cover * 100));
			printf("   Ice cover percentage:     %9.2f\n",
			       (node1->ice_cover * 100));
		}
		if (flag_moons && (node1->first_moon != NULL)) {
			printf("    MOONS:\n");
			printf("    #    Earth masses    orbital distance    radius    gravity\n");
			printf("                          (1000's of km)      (km)     (gees)\n");
			printf("    ---------------------------------------------------------------\n");
			for (node2=node1->first_moon, counter2=1;
			     node2 != NULL;
				node2=node2->next_planet, counter2++) {
				if (node2->mass_type == GAS_GIANT) {
					printf("    %2d   %2.2e            %5.2f        %3.1f     *gas giant*\n",
						  counter2,
						  (node2->mass*SUN_MASS_IN_EARTH_MASSES),
						  (node2->a * KM_PER_AU / 100000),
						  node2->radius);
				}
				else {
					printf("    %2d   %2.2e            %5.2f        %5.2f    %4.2f\n",
						  counter2,
						  (node2->mass*SUN_MASS_IN_EARTH_MASSES),
						  (node2->a * KM_PER_AU / 100000),
						  node2->radius,
						  node2->surf_grav);
				}
			}
		}
		else if (flag_moons && (node1->first_moon == NULL)) {
			printf("    NO MOONS\n");
		}
		printf("\n\n");
	}
}

void display_system(system)
sys_pointer system;
{
	if (flag_graphic)
		draw_system(system);
	else describe_system(system);
}

