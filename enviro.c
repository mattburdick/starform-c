/*----------------------------------------------------------------------*/
/*                              enviro.c                                */
/*                                                                      */
/* Planetary environment routines for determining things like planet    */
/* radius, density, mass, surface temp, etc.                            */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
#include    <math.h>
#include	<float.h>
#include    <errno.h>
#include	<stdio.h>

#include    "config.h"
#include    "const.h"
#include    "structs.h"
#include	"protos.h"

extern star_system sys;

extern int resonance;


/*--------------------------------------------------------------------------*/
/*   This function, given the orbital radius of a planet in AU, returns     */
/*   the orbital 'zone' of the particle.                                    */
/*--------------------------------------------------------------------------*/
int orb_zone(double orb_radius, double stell_luminosity_ratio)
{
    if (orb_radius < (4.0 * sqrt(stell_luminosity_ratio)))
        return(1);
    else {
        if ((orb_radius >= (4.0 * sqrt(stell_luminosity_ratio)))
            && (orb_radius < (15.0 * sqrt(stell_luminosity_ratio))))
            return(2);
        else return(3);
    }
}


/*--------------------------------------------------------------------------*/
/*   The mass is in units of solar masses, and the density is in units      */
/*   of grams/cc.  The radius returned is in units of km.                   */
/*--------------------------------------------------------------------------*/

double volume_radius(double mass, double density)
{
    double volume;

    mass = mass * SOLAR_MASS_IN_GRAMS;
    volume = mass / density;
    return(power((3.0 * volume) / (4.0 * PI),(1.0 / 3.0)) / CM_PER_KM);
}

/*--------------------------------------------------------------------------*/
/*   Returns the radius of the planet in kilometers.                        */
/*   The mass passed in is in units of solar masses.                        */
/*   This formula is listed as eq.9 in Fogg's article, although some typos  */
/*   crop up in that eq.  See "The Internal Constitution of Planets", by    */
/*   Dr. D. S. Kothari, Mon. Not. of the Royal Astronomical Society, vol 96 */
/*   pp.833-843, 1936 for the derivation.  Specifically, this is Kothari's  */
/*   eq.23, which appears on page 840.                                      */
/*--------------------------------------------------------------------------*/

double kothari_radius(double mass, int mass_type, int zone)
{
    double temp, temp2, atomic_weight, atomic_num;

    if (zone == 1) {
        if (mass_type == GAS_GIANT) {
            atomic_weight = 9.5;
            atomic_num = 4.5;
        }
        else {
            atomic_weight = 15.0;
            atomic_num = 8.0;
        }
    }
    else if (zone == 2) {
        if (mass_type == GAS_GIANT) {
            atomic_weight = 2.47;
            atomic_num = 2.0;
        }
        else {
            atomic_weight = 10.0;
            atomic_num = 5.0;
        }
    }
    else {
        if (mass_type == GAS_GIANT) {
            atomic_weight = 7.0;
            atomic_num = 4.0;
        }
        else {
            atomic_weight = 10.0;
            atomic_num = 5.0;
        }
    }
    temp = atomic_weight * atomic_num;
    temp = (2.0 * BETA_20 * power(SOLAR_MASS_IN_GRAMS,(1.0 / 3.0)))
        / (A1_20 * power(temp,(1.0 / 3.0)));
    temp2 = A2_20 * power(atomic_weight,(4.0 / 3.0))
        * power(SOLAR_MASS_IN_GRAMS,(2.0 / 3.0));
    temp2 = temp2 * power(mass,(2.0 / 3.0));
    temp2 = temp2 / (A1_20 * power(atomic_num, 2.0));
    temp2 = 1.0 + temp2;
    temp = temp / temp2;
    temp = (temp * power(mass,(1.0 / 3.0))) / CM_PER_KM;
    return(temp);
}


/*--------------------------------------------------------------------------*/
/*  The mass passed in is in units of solar masses, and the luminosity is   */
/*  is a unitless ratio.  The density is returned in units of grams/cc.     */
/*--------------------------------------------------------------------------*/

double empirical_density(double mass, double orb_radius, int mass_type, double star_luminosity)
{
    double temp, temp2;

    temp = power(mass * SUN_MASS_IN_EARTH_MASSES,(1.0 / 8.0));
    temp2 = sqrt(star_luminosity);
    temp = temp * power(temp2 / orb_radius, 0.25);
    if (mass_type == GAS_GIANT)
        return(temp * 1.2);
    else return(temp * 5.5);
}


/*--------------------------------------------------------------------------*/
/*  The mass passed in is in units of solar masses, and the equatorial      */
/*  radius is in km.  The density is returned in units of grams/cc.         */
/*--------------------------------------------------------------------------*/

double volume_density(double mass, double equat_radius)
{
    double volume;

    mass = mass * SOLAR_MASS_IN_GRAMS;
    equat_radius = equat_radius * CM_PER_KM;
    volume = (4.0 * PI * power(equat_radius, 3.0)) / 3.0;
    return(mass / volume);
}


/*--------------------------------------------------------------------------*/
/*  The separation is in units of AU, and both masses are in units of solar */
/*  masses.  The period returned is in terms of Earth days.                 */
/*--------------------------------------------------------------------------*/

double period(double separation, double small_mass, double large_mass)
{
    double period_in_years;

    period_in_years = sqrt(power(separation, 3.0)
                   / (small_mass + large_mass));
    return(period_in_years * DAYS_IN_A_YEAR);
}


/*--------------------------------------------------------------------------*/
/*   Fogg's information for this routine came from Dole "Habitable Planets  */
/* for Man", Blaisdell Publishing Company, NY, 1964.  From this, he came    */
/* up with his eq.12, which is the equation for the 'base_angular_velocity' */
/* below.  He then used an equation for the change in angular velocity per  */
/* time (dw/dt) from P. Goldreich and S. Soter's paper "Q in the Solar      */
/* System" in Icarus, vol 5, pp.375-389 (1966).  Using as a comparison the  */
/* change in angular velocity for the Earth, Fogg has come up with an       */
/* approximation for our new planet (his eq.13) and take that into account. */
/* This is used to find 'change_in_angular_velocity' below.                 */
/*                                                                          */
/*   Input parameters are mass (in solar masses), radius (in Km), orbital   */
/* period (in days), orbital radius (in AU), density (in g/cc),             */
/* eccentricity, and whether it is a gas giant or not.                      */
/*   The length of the day is returned in units of hours.                   */
/*--------------------------------------------------------------------------*/

double day_length(double mass, double radius, double eccentricity, double density,
          double orb_radius, double orb_period, int mass_type, double stell_mass_ratio, double age)
{
    double base_angular_velocity, planetary_mass_in_grams, k2,
    ang_velocity, equatorial_radius_in_cm, change_in_angular_velocity,
    spin_resonance_factor, year_in_hours, day_in_hours;
    int stopped = FALSE;

    resonance = FALSE;
    if (mass_type == GAS_GIANT)
        k2 = 0.24;
    else k2 = 0.33;
    planetary_mass_in_grams = mass * SOLAR_MASS_IN_GRAMS;
    equatorial_radius_in_cm = radius * CM_PER_KM;
    year_in_hours = orb_period * 24.0;
    base_angular_velocity = sqrt(2.0 * J * (planetary_mass_in_grams) /
                     (k2*power(equatorial_radius_in_cm,2.0)));
/*
 *  This next calculation determines how much the planet's rotation is
 *  slowed by the presence of the star.
 */
    change_in_angular_velocity = CHANGE_IN_EARTH_ANG_VEL
        * (density / EARTH_DENSITY)
        * (equatorial_radius_in_cm / EARTH_RADIUS)
        * (EARTH_MASS_IN_GRAMS / planetary_mass_in_grams)
        * power(stell_mass_ratio, 2.0)
        * (1.0 / power(orb_radius, 6.0));
    ang_velocity = base_angular_velocity + (change_in_angular_velocity * age);
/*
 *  Now we change from rad/sec to hours/rotation.
 */
    if (ang_velocity <= 0.0)
        stopped = TRUE;
    else day_in_hours = RADIANS_PER_ROTATION
        / (SECONDS_PER_HOUR * ang_velocity);
    if ((day_in_hours >= year_in_hours) || stopped) {
        if (eccentricity > 0.1) {
            spin_resonance_factor = (1.0 - eccentricity)
                                  / (1.0 + eccentricity);
            resonance = TRUE;
            return(spin_resonance_factor * year_in_hours);
        }
        else return(year_in_hours);
    }
    return(day_in_hours);
}


/*--------------------------------------------------------------------------*/
/*   The orbital radius is expected in units of Astronomical Units (AU).    */
/*   Inclination is returned in units of degrees.                           */
/*--------------------------------------------------------------------------*/

int inclination(double orb_radius)
{
    int temp;

    temp = (int)(power(orb_radius,0.2) * about(EARTH_AXIAL_TILT,0.4));
    return(temp % 360);
}


/*--------------------------------------------------------------------------*/
/*  This function implements the escape velocity calculation.  Note that    */
/*  it appears that Fogg's eq.15 is incorrect.                              */
/*  The mass is in units of solar mass, the radius in kilometers, and the   */
/*  velocity returned is in cm/sec.                                         */
/*--------------------------------------------------------------------------*/

double escape_vel(double mass, double radius)
{
    double mass_in_grams, radius_in_cm;

    mass_in_grams = mass * SOLAR_MASS_IN_GRAMS;
    radius_in_cm = radius * CM_PER_KM;
    return(sqrt(2.0 * GRAV_CONSTANT * mass_in_grams / radius_in_cm));
}


/*--------------------------------------------------------------------------*/
/*  This is Fogg's eq.16.  The molecular weight (usually assumed to be N2)  */
/*  is used as the basis of the Root Mean Square (RMS) velocity of the      */
/*  molecule or atom.  The velocity returned is in cm/sec.                  */
/*--------------------------------------------------------------------------*/

double rms_vel(double molecular_weight, double orb_radius, double luminosity)
{
    double exospheric_temp;

    exospheric_temp = EARTH_EXOSPHERE_TEMP
        * (luminosity / power(orb_radius, 2.0));
    return(sqrt((3.0 * MOLAR_GAS_CONST * exospheric_temp) / molecular_weight)
           * CM_PER_METER);
}


/*--------------------------------------------------------------------------*/
/*   This function returns the smallest molecular weight retained by the    */
/*  body, which is useful for determining the atmosphere composition.       */
/*  Orbital radius is in A.U.(ie: in units of the earth's orbital radius),  */
/*  mass is in units of solar masses, and equatorial radius is in units of  */
/*  kilometers.                                                             */
/*--------------------------------------------------------------------------*/

double molecule_limit(double mass, double equat_radius)
{
    double esc_velocity;

    esc_velocity = escape_vel(mass,equat_radius);
    return((3.0 * power(GAS_RETENTION_THRESHOLD * CM_PER_METER, 2.0)
        * MOLAR_GAS_CONST * EARTH_EXOSPHERE_TEMP)
        / power(esc_velocity,2.0));
}


/*--------------------------------------------------------------------------*/
/*   This function calculates the surface acceleration of a planet.  The    */
/*  mass is in units of solar masses, the radius in terms of km, and the    */
/*  acceleration is returned in units of cm/sec2.                           */
/*--------------------------------------------------------------------------*/

double accel(double mass, double radius)
{
    return(GRAV_CONSTANT * (mass * SOLAR_MASS_IN_GRAMS) /
           power(radius * CM_PER_KM, 2.0));
}


/*--------------------------------------------------------------------------*/
/*   This function calculates the surface gravity of a planet.  The         */
/*  acceleration is in units of cm/sec2, and the gravity is returned in     */
/*  units of Earth gravities.                                               */
/*--------------------------------------------------------------------------*/

double gravity(double acceleration)
{
    return(acceleration / EARTH_ACCELERATION);
}


/*--------------------------------------------------------------------------*/
/*  Note that if the orbital radius of the planet is greater than or equal  */
/*  to R_inner, 99% of it's volatiles are assumed to have been deposited in */
/*  surface reservoirs (otherwise, it suffers from the greenhouse effect).  */
/*--------------------------------------------------------------------------*/

int grnhouse(int zone, double orb_radius, double r_greenhouse)
{
    if ((orb_radius < r_greenhouse) && (zone == 1))
        return(TRUE);
    else return(FALSE);
}

/*--------------------------------------------------------------------------*/
/*  This implements Fogg's eq.17.  The 'inventory' returned is unitless.    */
/*  Returns a measure of the amount of gasses locked up inside the planet.  */
/*--------------------------------------------------------------------------*/

double vol_inventory(double mass, double esc_velocity, double rms_velocity, double stellar_mass, int zone,
             int greenhouse_effect)
{
    double velocity_ratio, proportion_const, temp1, temp2, earth_units;

    velocity_ratio = esc_velocity / rms_velocity;
    if (velocity_ratio >= GAS_RETENTION_THRESHOLD) {
        switch (zone) {
        case 1:
            proportion_const = 100000.0;
            break;
        case 2:
            proportion_const = 75000.0;
            break;
        case 3:
            proportion_const = 250.0;
            break;
        default:
            printf("ERROR: orb. zone not init'ed correctly!\n");
            break;
        }
        earth_units = mass * SUN_MASS_IN_EARTH_MASSES;
        temp1 = (proportion_const * earth_units) / stellar_mass;
        temp2 = about(temp1,0.2);
        if (greenhouse_effect)
            return(temp2);
        else return(temp2 / 100.0);
    }
    else return(0.0);
}


/*--------------------------------------------------------------------------*/
/*  This implements Fogg's eq.18, although it has been changed somewhat     */
/*  to account for planets so close to a star that their atmosphere has     */
/*  been blown off.                                                         */
/*  Input                   Units       Description                         */
/*  _____                   _____       ___________                         */
/*  volatile_gas_inventory  unitless    A measure of the gasses available   */
/*  equat_radius            kilometers  Equatorial radius of the planet     */
/*  grav                    gees        Surface gravity                     */
/*                                                                          */
/*  Output                  Units       Description                         */
/*  _____                   _____       ___________                         */
/*  pressure                millibars   Surface atmospheric pressure        */
/*--------------------------------------------------------------------------*/

double pressure(double volatile_gas_inventory, double equat_radius, double grav)
{
    equat_radius = KM_EARTH_RADIUS / equat_radius;
    return(volatile_gas_inventory * grav / power(equat_radius, 2.0));
}

/*--------------------------------------------------------------------------*/
/*   This function returns the boiling point of water in an atmosphere of   */
/*   pressure 'surf_pressure', given in millibars.  The boiling point is    */
/*   returned in units of Kelvin.  This is Fogg's eq.21.                    */
/*--------------------------------------------------------------------------*/

double boiling_point(double surf_pressure)
{
    double surface_pressure_in_bars;

    surface_pressure_in_bars = surf_pressure / MILLIBARS_PER_BAR;
    return(1.0 / (log(surface_pressure_in_bars) / -5050.5 + 1.0 / 373.0));
}


/*--------------------------------------------------------------------------*/
/*   This function is Fogg's eq.22.  Given the volatile gas inventory and   */
/*   planetary radius of a planet (in Km), this function returns the        */
/*   fraction of the planet covered with water.                             */
/*   I have changed the function very slightly:  the fraction of Earth's    */
/*   surface covered by water is 71%, not 75% as Fogg used.                 */
/*--------------------------------------------------------------------------*/

double hydro_fraction(double volatile_gas_inventory, double planet_radius)
{
    double temp;

    temp = (0.71 * volatile_gas_inventory / 1000.0)
        * power(KM_EARTH_RADIUS / planet_radius, 2.0);
    if (temp >= 1.0)
        return(1.0);
    else return(temp);
}


/*--------------------------------------------------------------------------*/
/*   Given the surface temperature of a planet (in Kelvin), this function   */
/*   returns the fraction of cloud cover available.  This is Fogg's eq.23.  */
/*   See Hart in "Icarus" (vol 33, pp23 - 39, 1978) for an explanation.     */
/*   This equation is Hart's eq.3.                                          */
/*   I have modified it slightly using constants and relationships from     */
/*   Glass's book "Introduction to Planetary Geology", p.46.                */
/*   The 'CLOUD_COVERAGE_FACTOR' is the amount of surface area on Earth     */
/*   covered by one Kg. of cloud.                                           */
/*--------------------------------------------------------------------------*/

double cloud_fraction(double surf_temp, double smallest_MW_retained, double equat_radius,
              double hyd_fraction)
{
    double water_vapor_in_kg, fraction, surf_area, hydro_mass;

    if (smallest_MW_retained > WATER_VAPOR)
        return(0.0);
    else {
        surf_area = 4.0 * PI * power(equat_radius, 2.0);
        hydro_mass = hyd_fraction * surf_area
            * EARTH_WATER_MASS_PER_AREA;
        water_vapor_in_kg = (0.00000001 * hydro_mass)
            * exp(Q2_36 * (surf_temp - 288.0));
        fraction = CLOUD_COVERAGE_FACTOR * water_vapor_in_kg / surf_area;
        if (fraction >= 1.0)
            return(1.0);
        else return(fraction);
    }
}


/*--------------------------------------------------------------------------*/
/*   Given the surface temperature of a planet (in Kelvin), this function   */
/*   returns the fraction of the planet's surface covered by ice.  This is  */
/*   Fogg's eq.24.  See Hart[24] in Icarus vol.33, p.28 for an explanation. */
/*   I have changed a constant from 70 to 90 in order to bring it more in   */
/*   line with the fraction of the Earth's surface covered with ice, which  */
/*   is approximatly .016 (=1.6%).                                          */
/*--------------------------------------------------------------------------*/

double ice_fraction(double hyd_fraction, double surf_temp)
{
    double temp;

    if (surf_temp > 328.0)
        surf_temp = 328.0;
    temp = power(((328.0 - surf_temp) / 90.0),5.0);
    if (temp > (1.5 * hyd_fraction))
        temp = (1.5 * hyd_fraction);
    if (temp >= 1.0)
        return(1.0);
    else return(temp);
}


/*--------------------------------------------------------------------------*/
/*  This is Fogg's eq.19.  The ecosphere radius is given in AU, the orbital */
/*  radius in AU, and the temperature returned is in Kelvin.                */
/*--------------------------------------------------------------------------*/

double eff_temp(double ecosphere_radius, double orb_radius, double albedo)
{
    return(sqrt(ecosphere_radius / orb_radius)
           * power((1.0 - albedo) / 0.7, 0.25)
           * EARTH_EFFECTIVE_TEMP);
}


/*--------------------------------------------------------------------------*/
/*  This is Fogg's eq.20, and is also Hart's eq.20 in his "Evolution of     */
/*  Earth's Atmosphere" article.  The effective temperature given is in     */
/*  units of Kelvin, as is the rise in temperature produced by the          */
/*  greenhouse effect, which is returned.                                   */
/*--------------------------------------------------------------------------*/

double green_rise(double optical_depth, double effective_temp, double surf_pressure)
{
    double convection_factor;

    convection_factor = EARTH_CONVECTION_FACTOR
        * power(surf_pressure / EARTH_SURF_PRES_IN_MILLIBARS, 0.25);
    return(power(1.0 + 0.75 * optical_depth, 0.25) - 1.0)
        * effective_temp * convection_factor;
}


/*--------------------------------------------------------------------------*/
/*   The surface temperature passed in is in units of Kelvin.               */
/*   The cloud adjustment is the fraction of cloud cover obscuring each     */
/*   of the three major components of albedo that lie below the clouds.     */
/*--------------------------------------------------------------------------*/

double planet_albedo(double water_fraction, double cld_fraction, double ice_frc, double surf_pressure)
{
    double rock_fraction, cloud_adjustment, components, cloud_part,
    rock_part, water_part, ice_part;

    rock_fraction = 1.0 - water_fraction - ice_frc;
    components = 0.0;
    if (water_fraction > 0.0)
        components = components + 1.0;
    if (ice_frc > 0.0)
        components = components + 1.0;
    if (rock_fraction > 0.0)
        components = components + 1.0;
    cloud_adjustment = cld_fraction / components;
    if (rock_fraction >= cloud_adjustment)
        rock_fraction = rock_fraction - cloud_adjustment;
    else rock_fraction = 0.0;
    if (water_fraction > cloud_adjustment)
        water_fraction = water_fraction - cloud_adjustment;
    else water_fraction = 0.0;
    if (ice_frc > cloud_adjustment)
        ice_frc = ice_frc - cloud_adjustment;
    else ice_frc = 0.0;
    cloud_part = cld_fraction * about(CLOUD_ALBEDO,0.2);
    if (surf_pressure == 0.0)
        rock_part = rock_fraction * about(ROCKY_AIRLESS_ALBEDO,0.3);
    else rock_part = rock_fraction * about(ROCKY_ALBEDO,0.1);
    water_part = water_fraction * about(WATER_ALBEDO,0.2);
    if (surf_pressure == 0.0)
        ice_part = ice_frc * about(AIRLESS_ICE_ALBEDO,0.4);
    else ice_part = ice_frc * about(ICE_ALBEDO,0.1);
    return(cloud_part + rock_part + water_part + ice_part);
}


/*--------------------------------------------------------------------------*/
/*   This function returns the dimensionless quantity of optical depth,     */
/*   which is useful in determining the amount of greenhouse effect on a    */
/*   planet.                                                                */
/*--------------------------------------------------------------------------*/

double opacity(double molecular_weight, double surf_pressure)
{
    double optical_depth;

    optical_depth = 0.0;
    if ((molecular_weight >= 0.0) && (molecular_weight < 10.0))
        optical_depth = optical_depth + 3.0;
    if ((molecular_weight >= 10.0) && (molecular_weight < 20.0))
        optical_depth = optical_depth + 2.34;
    if ((molecular_weight >= 20.0) && (molecular_weight < 30.0))
        optical_depth = optical_depth + 1.0;
    if ((molecular_weight >= 30.0) && (molecular_weight < 45.0))
        optical_depth = optical_depth + 0.15;
    if ((molecular_weight >= 45.0) && (molecular_weight < 100.0))
        optical_depth = optical_depth + 0.05;
    if (surf_pressure >= (70.0 * EARTH_SURF_PRES_IN_MILLIBARS))
        optical_depth = optical_depth * 8.333;
    else if (surf_pressure >= (50.0 * EARTH_SURF_PRES_IN_MILLIBARS))
        optical_depth = optical_depth * 6.666;
    else if (surf_pressure >= (30.0 * EARTH_SURF_PRES_IN_MILLIBARS))
        optical_depth = optical_depth * 3.333;
    else if (surf_pressure >= (10.0 * EARTH_SURF_PRES_IN_MILLIBARS))
        optical_depth = optical_depth * 2.0;
    else if (surf_pressure >= (5.0 * EARTH_SURF_PRES_IN_MILLIBARS))
        optical_depth = optical_depth * 1.5;
    return(optical_depth);
}


/*--------------------------------------------------------------------------*/
/*   The temperature calculated is in degrees Kelvin.                       */
/*   Quantities already known which are used in these calculations:         */
/*     planet->molec_weight                                                 */
/*     planet->surf_pressure                                                */
/*       R_ecosphere                                                        */
/*     planet->a                                                            */
/*     planet->volatile_gas_inventory                                       */
/*     planet->radius                                                       */
/*     planet->boil_point                                                   */
/*  The 'counter' variable used in the interation loop is used to break     */
/*  out of the loop after 100 iterations - just in case the temperature     */
/*  refuses to converge.                                                    */
/*--------------------------------------------------------------------------*/

void iterate_surface_temp(planet_pointer *planet, double r_ecosphere)
{
    double effective_temp, greenhs_rise, previous_temp,
    optical_depth, albedo, water, eff_water, clouds, ice, new_temp = 0.0;
    int counter = 0;

    albedo = EARTH_ALBEDO;
    water = hydro_fraction((*planet)->volatile_gas_inventory,
                   (*planet)->radius);
    optical_depth = opacity((*planet)->molec_weight,
                (*planet)->surf_pressure);
    do {
        effective_temp = eff_temp(r_ecosphere,
                      (*planet)->a,
                      albedo);
        if (counter == 0) {
            previous_temp = effective_temp;
        } else {
            previous_temp = new_temp;
        }
        greenhs_rise = green_rise(optical_depth,
                      effective_temp,
                      (*planet)->surf_pressure);
		/*
		 *	Watch out for overflow:
		 */
		if (effective_temp + greenhs_rise > DBL_MAX) {
			effective_temp = DBL_MAX;
			greenhs_rise = 0.0;
			new_temp= DBL_MAX;
		}
		else new_temp = effective_temp + greenhs_rise;
        clouds = cloud_fraction(new_temp,
                    (*planet)->molec_weight,
                    (*planet)->radius,
                    water);
        ice = ice_fraction(water,new_temp);
        if ((new_temp >= (*planet)->boil_point)
            || (new_temp <= FREEZING_POINT_OF_WATER))
            eff_water = 0.0;
        else eff_water = water;
        albedo = planet_albedo(eff_water,
                       clouds,
                       ice,
                       (*planet)->surf_pressure);
        counter++;
    } while ((fabs(new_temp - previous_temp) > 1.0)
          && (counter < TEMP_ITERATION_LIMIT));
    (*planet)->hydrosphere = eff_water;

#ifdef NOTDEF
    optical_depth = opacity((*planet)->molec_weight,
                (*planet)->surf_pressure);
    effective_temp = eff_temp(r_ecosphere,(*planet)->a,EARTH_ALBEDO);
    greenhs_rise = green_rise(optical_depth,
                  effective_temp,
                  (*planet)->surf_pressure);
    new_temp = effective_temp + greenhs_rise;
    previous_temp = new_temp - 5.0;
    while ((fabs(new_temp - previous_temp) > 1.0) && (counter < TEMP_ITERATION_LIMIT)) {
        counter++;
        previous_temp = new_temp;
        water = hydro_fraction((*planet)->volatile_gas_inventory,
                       (*planet)->radius);
        clouds = cloud_fraction(new_temp,
                    (*planet)->molec_weight,
                    (*planet)->radius,
                    water);
        ice = ice_fraction(water,new_temp);
        if ((new_temp >= (*planet)->boil_point)
            || (new_temp <= FREEZING_POINT_OF_WATER))
            water = 0.0;
        albedo = planet_albedo(water,
                       clouds,
                       ice,
                       (*planet)->surf_pressure);
        optical_depth = opacity((*planet)->molec_weight,
                    (*planet)->surf_pressure);
        effective_temp = eff_temp(r_ecosphere,(*planet)->a,albedo);
        greenhs_rise = green_rise(optical_depth,
                      effective_temp,
                      (*planet)->surf_pressure);
        new_temp = effective_temp + greenhs_rise;
    }
    (*planet)->hydrosphere = water;
#endif
    (*planet)->cloud_cover = clouds;
    (*planet)->ice_cover = ice;
    (*planet)->albedo = albedo;
    (*planet)->surf_temp = new_temp;
}

