/*----------------------------------------------------------------------*/
/*                               protos.h                               */
/*                                                                      */
/*  Header file containing all function prototypes.                     */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
#include "structs.h"
/*
 *	From 'accrete.c':
 */
void           set_initial_conditions(double, double);
planet_pointer sorted_list_insert(planet_pointer, planet_pointer);
double         stell_dust_limit(double, double, int);
double         nearest_body(double);
double         farthest_body(double);
double         roche_limit(double);
double         inner_effect_limit(double, double, double);
double         outer_effect_limit(double, double, double);
int            dust_available(double, double);
dust_pointer   prior_dust_band(dust_pointer, dust_pointer);
double         collect_dust(double, double, double, double, dust_pointer);
double         critical_limit(double, double, double);
double         accrete_dust(double, double, double, double);
planet_pointer find_collision (planet_pointer, double, double);
void           collide_planets(double, double, double, planet_pointer, double);
void           coalesce_planetesimals(double, double, double, double, double, int);
planet_pointer dist_masses(double, double, int, planet_pointer, double);
planet_pointer check_planets(planet_pointer, double, double);
planet_pointer init_planet_list(star_pointer);

/*
 *	From 'utils.c':
 */
double         power(double, double);
double         random_number(double, double);
double         about(double, double);
double         random_eccentricity(void);

/*
 *	From 'enviro.c':
 */
int            orb_zone(double, double);
double         volume_radius(double, double);
double         kothari_radius(double, int, int);
double         empirical_density(double, double, int, double);
double         volume_density(double, double);
double         period(double, double, double);
double         day_length(double, double, double, double, double, double, int, double, double);
int            inclination(double);
double         escape_vel(double, double);
double         rms_vel(double, double, double);
double         molecule_limit(double, double);
double         accel(double, double);
double         gravity(double);
int            grnhouse(int, double, double);
double         vol_inventory(double, double, double, double, int, int);
double         pressure(double, double, double);
double         boiling_point(double);
double         hydro_fraction(double, double);
double         cloud_fraction(double, double, double, double);
double         ice_fraction(double, double);
double         eff_temp(double, double, double);
double         green_rise(double, double, double);
double         planet_albedo(double, double, double, double);
double         opacity(double, double);
void           iterate_surface_temp(planet_pointer *, double);

/*
 *	From 'starform.c':
 */
void           usage(char *);
void           init(void);
void           generate_stellar_system(void);
void           main(int, char **);

/*
 *	From 'stars.c':
 */
double         luminosity(double, int);
double         star_radius(double, int, int);
double         star_age(double);
char *         classify(double, int);
double         star_mass(int, char, int);
int            verify_startype(char, int, char);
void           startype_error(int, char, int, char);
double         rand_star_mass(int);
int            rand_type(void);

/*
 *	From 'display.c':
 */
void           draw_system(sys_pointer);
void           describe_star(star_pointer);
void           describe_system(sys_pointer);
void           display_system(sys_pointer);

