/*----------------------------------------------------------------------*/
/*                              accrete.c                               */
/*                                                                      */
/* Planetary accretion routines based on:                               */
/*                                                                      */
/* Dole, Stephen H.  "Formation of Planetary Systems by Aggregation:    */
/* a Computer Simulation"  October 1969,  Rand Corporation Paper P-4226 */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/

#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#include    <errno.h>

#include    "config.h"
#include    "const.h"
#include    "structs.h"
#include	"protos.h"

extern int flag_verbose;

/*
 * A few variables global to the entire program:
 */
planet_pointer planet_head;

/*
 * Now for some variables global to the accretion process:
 */
int dust_left;
double r_inner, r_outer, reduced_mass, dust_density;
dust_pointer dust_head;


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void set_initial_conditions(double inner_limit_of_dust, double outer_limit_of_dust)
{
    if ((dust_head = (dust *)malloc((unsigned)sizeof(dust))) == NULL) {
        perror("malloc'ing head of dust list");
        exit(1);
    }
    dust_head->next_band = NULL;
    dust_head->outer_edge = outer_limit_of_dust;
    dust_head->inner_edge = inner_limit_of_dust;
    dust_head->dust_present = TRUE;
    dust_head->gas_present = TRUE;
    dust_left = TRUE;
	if (flag_verbose >= LEVEL3) {
		printf("      Creating the head of the dust list (%4.2lg - %4.2lg).\n",
			dust_head->inner_edge, dust_head->outer_edge);
	}
}

/*--------------------------------------------------------------------------*/
/*  Insert the given planet into a list of planets sorted by distance from  */
/*  the primary.                                                            */
/*--------------------------------------------------------------------------*/
planet_pointer sorted_list_insert (planet_pointer head, planet_pointer planet)
{
    planet_pointer node;
    planet_pointer trailer = NULL;

    /*
     *  If the list was empty to begin with, just return the planet:
     */
    if (head == NULL) {
        planet->next_planet = NULL;
        return(planet);
    }
    /*
     *  Traverse the list until we find one further out than 'planet',
     *  then insert 'planet' in front of it:
     */
    for (node=head; (node != NULL); trailer=node, node=node->next_planet) {
        if (node->a > planet->a) {
            if (trailer == NULL) {
                planet->next_planet = node;
                return(planet);
            }
            else {
                planet->next_planet = trailer->next_planet;
                trailer->next_planet = planet;
                return(head);
            }
        }
    }
    /*
     *  If we've gotten this far, the planet must be inserted at the
     *  end of the list:
     */
    trailer->next_planet = planet;
    planet->next_planet = NULL;
    return(head);
}

/*
 *  primary_effect is the term in the equation that reduces the size of the
 *  dust limit due to the proximity of the primary star.
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double stell_dust_limit(double mass_ratio, double dist_from_primary, int central_mass)
{
    double norm_limit, primary_effect;

    norm_limit= 200.0 * power(mass_ratio,(1.0 / 3.0));
    if (central_mass == PLANET) {
        norm_limit = norm_limit / 125.0;
        primary_effect = power(dist_from_primary, 2.0);
        if (primary_effect <= 1.0)
            return (norm_limit * primary_effect);
        else return (norm_limit);
    }
    else return(norm_limit);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double nearest_body(double mass_ratio)
{
    return(0.3 * power(mass_ratio,(1.0 / 3.0)));
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double farthest_body(double stell_mass_ratio)
{
    return(50.0 * power(stell_mass_ratio,(1.0 / 3.0)));
}

/*--------------------------------------------------------------------------*/
/*   This function returns the distance from a planet at which a moon can   */
/*   no longer hold itself together through its gravitational field.        */
/*   Inside this limit, the moon will eventually break up from tidal action.*/
/*   The input diameter is in units of Km, so we must first convert to AU   */
/*   then multiply it by Roche's limit.  The output is in units of AU.      */
/*--------------------------------------------------------------------------*/
double roche_limit(double diameter)
{
    double dia_in_AU;

    dia_in_AU = diameter / KM_PER_AU;
    return(2.44 * dia_in_AU);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double inner_effect_limit(double a, double e, double mass)
{
    return (a * (1.0 - e) * (1.0 - mass) / (1.0 + CLOUD_ECCENTRICITY));
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
double outer_effect_limit(double a, double e, double mass)
{
    return (a * (1.0 + e) * (1.0 + mass) / (1.0 - CLOUD_ECCENTRICITY));
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
int dust_available(double inside_range, double outside_range)
{
    dust_pointer current_dust_band;

    current_dust_band = dust_head;
    while ((current_dust_band != NULL)
           && (current_dust_band->outer_edge < inside_range)) {
        current_dust_band = current_dust_band->next_band;
    }
    if (current_dust_band == NULL) {
		return(FALSE);
    }
    else if (current_dust_band->dust_present) {
    	return(TRUE);
    }
    while ((current_dust_band != NULL)
           && (current_dust_band->inner_edge < outside_range)) {
		if (current_dust_band->dust_present) {
			return(TRUE);
		}
        current_dust_band = current_dust_band->next_band;
    }
    return(FALSE);
}

/*--------------------------------------------------------------------------*/
/*  Given a pointer to a dust band, this function returns the dust band     */
/*  immediately preceding it in the list.  If it was the first in the list, */
/*  a NULL will be returned.                                                */
/*--------------------------------------------------------------------------*/
dust_pointer prior_dust_band (dust_pointer dust_head, dust_pointer band)
{
    dust_pointer temp, trailer;

    trailer = NULL;
    for (temp = dust_head; (temp != NULL); trailer = temp, temp = temp->next_band) {
        if (temp == band) {
            return(trailer);
        }
    }
    printf("ERROR: searching for nonexistant dust band\n");
    exit(1);
}

/*--------------------------------------------------------------------------*/
/*  This routine compares the location of a test mass with the location     */
/*  of any dust and gas bands remaining.  Any dust band that lies within    */
/*  the object's range of gravitational effect (from 'r_inner' to 'r_outer')*/
/*  is swept up by the object.  Additionally, if the object's mass is       */
/*  greater than the critical mass, gas is also collected by the object.    */
/*                                                                          */
/*  The new mass for the object is returned.                                */
/*                                                                          */
/*  VARIABLES PASSED IN:                                                    */
/*    mass         Mass of the accreting object (Solar masses)              */
/*    a            Distance from the primary star (AUs)                     */
/*    e            Eccentricity of the object's orbit                       */
/*    crit_mass    Mass at which, for this orbit and star, a normal planet  */
/*                 begins to sweep up gas as well as dust and become a      */
/*                 gas giant.                                               */
/*    dust_head    Pointer to the head of the dust band list                */
/*                                                                          */
/*  LOCAL VARIABLES:                                                        */
/*    r_inner      Innermost gravitational effect limit of the object       */
/*    r_outer      Outermost gravitational effect limit of the object       */
/*                                                                          */
/*--------------------------------------------------------------------------*/
double collect_dust(double mass, double a, double e, double crit_mass, dust_pointer dust_head)
{
    double mass_density, temp1, temp2, bandwidth, width, volume,
            accumulated_mass;
    dust_pointer band, newband, gasband, prev_band;

    /*
     *  Find the effective mass and its range of effect ('r_inner'
     *  through 'r_outer'):
     */
    reduced_mass = power(mass / (1.0 + mass), 0.25);
    r_inner = inner_effect_limit(a, e, reduced_mass);
    r_outer = outer_effect_limit(a, e, reduced_mass);
    if (r_inner < 0.0) {
        printf("error!\n");
		exit(1);
    }
    /*
     *  Visit each dust band and add up any dust collected from each.
     *  Start with the original mass of the object:
     */
    accumulated_mass = mass;
    for (band = dust_head; (band != NULL); band = band->next_band) {
        /*
         *  If there is no gas in this band OR if the band lies outside
         *  the range of effect completely OR if no dust is present and
         *  the mass is too small to pick up gas, go to the next band:
         */
        if (band->gas_present == FALSE) {
            continue;
        }
        if ((band->outer_edge <= r_inner) || (band->inner_edge >= r_outer)) {
            continue;
        }
        if ((mass < crit_mass) && (band->dust_present == FALSE)) {
            continue;
        }
        /*
         *  Dust or gas exists and lies within range - sweep some up:
         */
        if (mass < crit_mass) {
            mass_density = dust_density;
        }
        else {
            mass_density = K * dust_density /
                (1.0 + sqrt(crit_mass / mass) * (K - 1.0));
        }
        bandwidth = (r_outer - r_inner);
        temp1 = r_outer - band->outer_edge;
        if (temp1 < 0.0)
            temp1 = 0.0;
        temp2 = band->inner_edge - r_inner;
        if (temp2 < 0.0)
            temp2 = 0.0;
        width = bandwidth - temp1 - temp2;
        volume = 4.0 * PI * power(a,2.0) * reduced_mass
            * (1.0 - e * (temp1 - temp2) / bandwidth) * width;
        accumulated_mass = accumulated_mass + (volume * mass_density);
    }
	/*
	 *	Now re-visit each band removing dust and gas if necessary.  It
	 *	may also be necessary to reduce the size of some dust bands
	 *	and create new gas bands (if the dust is removed from a dust
	 *	band, it becomes a gas band).
	 */
	for (band = dust_head; (band != NULL); band = band->next_band) {
        if (band->gas_present == FALSE) {
            continue;
        }
        if ((band->outer_edge <= r_inner) || (band->inner_edge >= r_outer)) {
            continue;
        }
        if ((mass < crit_mass) && (band->dust_present == FALSE)) {
            continue;
        }
        temp1 = r_outer - band->outer_edge;
        if (temp1 < 0.0)
            temp1 = 0.0;
        temp2 = band->inner_edge - r_inner;
        if (temp2 < 0.0)
            temp2 = 0.0;
        /*
         *  Some dust has been swept up, so update this band:
         */
        /*
         *  Case 1: the area of effect lies entirely within the dust band:
         *  Result: divide the original dust band into two smaller ones:
         */
        if ((temp1 == 0.0) && (temp2 == 0.0)) {
            if ((newband = (dust *)malloc((unsigned)sizeof(dust))) == NULL) {
                perror("creating a new dust band");
                exit(1);
            }
            newband->inner_edge = r_outer;
            newband->outer_edge = band->outer_edge;
            newband->dust_present = band->dust_present;
            newband->gas_present = band->gas_present;
            newband->next_band = band->next_band;
            band->outer_edge = r_inner;
            band->next_band = newband;
            if (flag_verbose >= LEVEL3) {
                printf("      Creating a new dust band 1 (%4.2lg - %4.2lg).\n",
            		newband->inner_edge, newband->outer_edge);
            }
            if (mass < crit_mass) {
                /*
                 *  The mass isn't a gas giant, so it'll sweep away all the
                 *  dust in it's range, but leave the gas.  Therefore, we
                 *  need to create a new gas band here:
                 */
                if ((gasband = (dust *)malloc((unsigned)sizeof(dust))) == NULL) {
                    perror("creating a new gas band");
                    exit(1);
                }
                gasband->inner_edge = r_inner;
                gasband->outer_edge = r_outer;
                gasband->dust_present = FALSE;
                gasband->gas_present = TRUE;
                gasband->next_band = newband;
                band->next_band = gasband;
                if (flag_verbose >= LEVEL3) {
                    printf("      Creating a new gas band 2 (%4.2lg - %4.2lg).\n",
            		gasband->inner_edge, gasband->outer_edge);
                }
            }
        }
        /*
         *  Case 2: the area of effect encompasses the dust band entirely:
         *  Result: Remove the band if both dust and gas can be swept up,
         *          otherwise, just remove the dust from it:
         */
        else if ((temp1 > 0.0) && (temp2 > 0.0)) {
            if (mass >= crit_mass) {
                prev_band = prior_dust_band(dust_head, band);
                if (prev_band == NULL) {
                    dust_head = band->next_band;
                }
                else {
                    prev_band->next_band = band->next_band;
                }
                if (flag_verbose >= LEVEL3) {
                    printf("      Freeing a gas band 3 (%4.2lg - %4.2lg).\n",
	                	band->inner_edge, band->outer_edge);
                }
                free(band);
            }
            else {
                band->dust_present = FALSE;
                if (flag_verbose >= LEVEL3) {
                    printf("      Removing dust from a dust/gas band 4 (%4.2lg - %4.2lg).\n",
	                	band->inner_edge, band->outer_edge);
                }
            }
        }
        /*
         *  Case 3: the area of effect and the dust band overlap with
         *          the dust band slightly further from the primary star.
         *  Result: Remove the inner part of the band if both dust and gas
         *          can be swept up, otherwise, just remove the dust from it:
         */
        else if (temp2 > 0.0) {
            if (mass >= crit_mass) {
                band->inner_edge = r_outer;
                if (flag_verbose >= LEVEL3) {
                    printf("      Reducing a gas band 5 (%4.2lg - %4.2lg).\n",
	                	band->inner_edge, band->outer_edge);
                }
            }
            else {
                /*
                 *  If there is a dust band prior to this one in the list
                 *  and it only has gas in it already and touches the inner
                 *  edge of the current band, we don't need to create a
                 *  new band - just add the current one onto the prior one.
                 */
                prev_band = prior_dust_band(dust_head, band);
                if (prev_band != NULL) {
                    if ((prev_band->dust_present == FALSE)
                        && (prev_band->outer_edge == band->inner_edge)) {
                        prev_band->outer_edge = r_outer;
                        band->inner_edge = r_outer;
                        if (flag_verbose >= LEVEL3) {
                        	printf("      Increasing a gas band 6 (%4.2lg - %4.2lg).\n",
	                        	prev_band->inner_edge, prev_band->outer_edge);
                        	printf("      Reducing a dust band 7 (%4.2lg - %4.2lg).\n",
	                        	band->inner_edge, band->outer_edge);
                        }
                        continue;
                    }
                }
                if ((gasband = (dust *)malloc((unsigned)sizeof(dust))) == NULL) {
                    perror("creating a new gas band");
                    exit(1);
                }
                gasband->inner_edge = band->inner_edge;
                gasband->outer_edge = r_outer;
                gasband->dust_present = FALSE;
                gasband->gas_present = TRUE;
                gasband->next_band = band;
                band->inner_edge = r_outer;
                if (flag_verbose >= LEVEL3) {
                	printf("      Reducing a dust band 8 (%4.2lg - %4.2lg).\n",
                		band->inner_edge, band->outer_edge);
                    printf("      Creating a new gas band 9 (%4.2lg - %4.2lg).\n",
	                	gasband->inner_edge, gasband->outer_edge);
                }
                prev_band = prior_dust_band(dust_head, band);
                if (prev_band == NULL) {
                    dust_head = gasband;
                }
                else {
                    prev_band->next_band = gasband;
                }
            }
        }
        /*
         *  Case 4: the area of effect and the dust band overlap with
         *          the dust band slightly closer to the primary star.
         *  Result: Remove the outer part of the band if both dust and gas
         *          can be swept up, otherwise, just remove the dust from it:
         */
        else if (temp1 > 0.0) {
            if (mass >= crit_mass) {
                band->outer_edge = r_inner;
            }
            else {
                /*
                 *  As above, if there is a dust band after this one
                 *  and it only has gas in it already and touches the outer
                 *  edge of the current band, we don't need to create a
                 *  new band - just add the current one onto the next one.
                 */
                if (band->next_band != NULL) {
                    if ((band->next_band->dust_present == FALSE)
                        && (band->next_band->inner_edge == band->outer_edge)) {
                        band->next_band->inner_edge = r_inner;
                        band->outer_edge = r_inner;
                        if (flag_verbose >= LEVEL3) {
                        	printf("      Increasing a gas band 10 (%4.2lg - %4.2lg).\n",
	                        	band->next_band->inner_edge, band->next_band->outer_edge);
                        	printf("      Reducing a dust band 11 (%4.2lg - %4.2lg).\n",
	                        	band->inner_edge, band->outer_edge);
                        }
                        continue;
                    }
                }
                if ((gasband = (dust *)malloc((unsigned)sizeof(dust))) == NULL) {
                    perror("creating a new gas band");
                    exit(1);
                }
                gasband->inner_edge = r_inner;
                gasband->outer_edge = band->outer_edge;
                gasband->dust_present = FALSE;
                gasband->gas_present = TRUE;
                gasband->next_band = band->next_band;
                band->outer_edge = r_inner;
                band->next_band = gasband;
                if (flag_verbose >= LEVEL3) {
                	printf("      Reducing a dust band 12 (%4.2lg - %4.2lg).\n",
                		band->inner_edge, band->outer_edge);
                    printf("      Creating a new gas band 13 (%4.2lg - %4.2lg).\n",
	                	gasband->inner_edge, gasband->outer_edge);
                }
            }
        }
    }
    return(accumulated_mass);
}


/*--------------------------------------------------------------------------*/
/*   Orbital radius is in AU, eccentricity is unitless, and the stellar     */
/*  luminosity ratio is with respect to the sun.  The value returned is the */
/*  mass at which the planet begins to accrete gas as well as dust, and is  */
/*  in units of solar masses.                                               */
/*--------------------------------------------------------------------------*/

double critical_limit(double orb_radius, double eccentricity, double stell_luminosity_ratio)
{
    double temp, perihelion_dist;

    perihelion_dist = (orb_radius - orb_radius * eccentricity);
    temp = perihelion_dist * sqrt(stell_luminosity_ratio);
    return(B * power(temp,-0.75));
}



/*--------------------------------------------------------------------------*/
/*  Given a mass at a particular orbit, this function repeatedly calls      */
/*  'collect_dust' to sweep up any dust and gas it can.  Each successive    */
/*  call to 'collect_dust' is done with the original mass plus additional   */
/*  mass from sweeping up dust previously.  The process stops when the mass */
/*  accumulation slows.                                                     */
/*--------------------------------------------------------------------------*/
double accrete_dust(double mass, double a, double e, double crit_mass)
{
    double new_mass;
	dust_pointer band;

    new_mass = mass;
    do {
        mass = new_mass;
        new_mass = collect_dust(new_mass, a, e, crit_mass, dust_head);
    }
    while ((new_mass - mass) > (0.001 * mass));
	/*
	 *  Traverse the dust bands to check if there is any dust remaining.
	 *  The global boolean 'dust_left' is used in 'dist_masses'.
	 */
	dust_left = FALSE;
    for (band = dust_head; (band != NULL); band = band->next_band) {
        if (band->dust_present) {
            dust_left = TRUE;
			break;
        }
    }
    return(new_mass);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
planet_pointer find_collision (planet_pointer head, double a, double e)
{
    planet_pointer node;
    planet_pointer closest_neighbor = NULL;
    double         closest_approach = 0.0;
    double         separation, dist1, dist2;

    for (node = head; (node); node = node->next_planet) {
        separation = node->a - a;
        /*
         *  In the following calculations, 'dist1' is the distance over
         *  which the new planet gravitationally attracts the existing
         *  planet while 'dist2' is the distance over which the existing
         *  planet affects the new one.  A collision occurs if the
         *  separation of the two planets is less than the gravitational
         *  effects distance of either.
         */
        reduced_mass = power(node->mass / (1.0 + node->mass),0.25);
        if ((separation > 0.0)) {
            /*
             *  The neighbor is farther from the star than our test planet:
             */
            dist1 = (a * (1.0 + e) * (1.0 + reduced_mass)) - a;
            dist2 = node->a - (node->a*(1.0 - node->e)*(1.0 - reduced_mass));
        }
        else {
            /*
             *  The new planet is farther from the star than it's neighbor:
             */
            dist1 = a - (a * (1.0 - e) * (1.0 - reduced_mass));
            dist2 = (node->a*(1.0 + node->e)*(1.0 + reduced_mass)) - node->a;
        }

        if ((   (fabs(separation) <= fabs(dist1))
             || (fabs(separation) <= fabs(dist2)))) {
            /*
             *  These two should collide.  Save the distance and check the
             *  next planet.  If it is farther away, stop looking for
             *  other collisions and return the saved one:
             */
            if (closest_neighbor) {
                /*
                 *  Already found a planet to collide with - is this
                 *  one closer?
                 */
                if (fabs(separation) < closest_approach) {
                    closest_approach = fabs(separation);
                    closest_neighbor = node;
                }
            }
            else {
                /*
                 *  No other planet is close enough so far, so save the
                 *  current info and continue checking:
                 */
                closest_neighbor = node;
                closest_approach = fabs(separation);
            }
        }
    }
    return(closest_neighbor);
}


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
void collide_planets(double a, double e, double mass, planet_pointer node, double stell_luminosity_ratio)
{
    double          new_orbit, temp;

    new_orbit = (node->mass + mass) /((node->mass / node->a) + (mass / a));

    if (flag_verbose >= LEVEL1) {
        switch (node->mass_type) {
        case STAR:
            printf("  Collision with a star! (%4.2lg, %4.2lg -> %4.2lg)\n",
                   a, node->a, new_orbit);
            break;
        case PLANET:
            printf("  Collision with a planet! (%4.2lg, %4.2lg -> %4.2lg)\n",
                   a, node->a, new_orbit);
            break;
        case GAS_GIANT:
            printf("  Collision with a gas giant! (%4.2lg, %4.2lg -> %4.2lg)\n",
                   a, node->a, new_orbit);
            break;
        case MOON:
            printf("  Collision with a moon! (%4.2lg, %4.2lg -> %4.2lg)\n",
                   a, node->a, new_orbit);
            break;
        default:
            printf("ERROR: collision with unknown mass type %d\n",
                node->mass_type);
            exit(1);
            break;
        }
    }

    temp = node->mass * sqrt(node->a) * sqrt(1.0 - power(node->e,2.0));
    temp = temp + (mass * sqrt(a) * sqrt(sqrt(1.0 - power(e,2.0))));
    temp = temp / ((node->mass + mass) * sqrt(new_orbit));
    temp = 1.0 - power(temp,2.0);
    if (((temp < 0.0) || (temp >= 1.0))) {
        temp = 0.0;
    }

    e = sqrt(temp);
    temp = node->mass + mass;
    temp = accrete_dust(temp,new_orbit,e,stell_luminosity_ratio);

    node->a = new_orbit;
    node->e = e;
    node->mass = temp;

    /*
     *  If the protoplanet had the misfortune to collide with a
     *  star, update the corresponding star node for the node
     *  planet node:
     */
    if (node->mass_type == STAR) {
        node->star_ptr->orbit_radius = node->a;
        node->star_ptr->stell_mass_ratio = node->mass;
    }
}


/*--------------------------------------------------------------------------*/
/*  coalesce_planetesimals checks if the protoplanet described by a, e,     */
/*  mass, etc, crosses the orbits of any planets already generated.  If so, */
/*  the masses of the two planets are added (therefore assuming a perfectly */
/*  inelastic collision) and an orbit for the resulting new planet          */
/*  computed.  The new mass is then allowed to accrete more dust.           */
/*  If no collision with another planet occurrs, a new planet is created    */
/*  its statistics filled in with those of the protoplanet.                 */
/*--------------------------------------------------------------------------*/
void coalesce_planetesimals(double a, double e, double mass, double crit_mass,
                double stell_luminosity_ratio, int orbit_type)
{
    planet_pointer node, new_planet;
    int finished;

    finished = FALSE;
    if (mass <= TRIVIAL_MASS) {
        if (flag_verbose >= LEVEL1) {
            printf("  Trivial mass (%9.3f Earth masses) - not adding it.\n",
        			mass * SUN_MASS_IN_EARTH_MASSES);
        }
        return;
    }
    node = find_collision(planet_head, a, e);
    if (node) {
        /*
         *  node is the closest planet in a colliding orbit with the
         *  planet described by a, e, and mass.
         */
        collide_planets(a, e, mass, node, stell_luminosity_ratio);
    }
    else {
        /*
         *  The new planet won't collide with any other planet or star,
         *  so allocate space for it and insert it into the system's
         *  linked list:
         */
        if ((new_planet = (planets *)malloc((unsigned)sizeof(planets)))
            == NULL) {
            perror("malloc'ing a new planet");
            exit(1);
        }
        if (flag_verbose >= LEVEL3) {
            printf("      Creating a new planet.\n");
        }
        new_planet->mass_type = orbit_type;
        new_planet->a = a;
        new_planet->e = e;
        new_planet->first_moon = NULL;
        new_planet->next_planet = NULL;
        if ((mass >= crit_mass)) {
            new_planet->mass_type = GAS_GIANT;
        }
        else {
            new_planet->mass_type = orbit_type;
        }
        new_planet->mass = mass;

        planet_head = sorted_list_insert(planet_head, new_planet);
    }
}


/*--------------------------------------------------------------------------*/
/*  This function builds a linked list of planets by repeatedly injecting   */
/*  protoplanets into the gas and dust cloud about a star.  The process     */
/*  will end when all of the dust has been swept up by the planets.         */
/*  'orbit_type' may be either PLANET (indicating that we're building a     */
/*  series of planetary bodies about a star) or MOON (indicating we're      */
/*  building moons around a planet).                                        */
/*--------------------------------------------------------------------------*/
planet_pointer dist_masses(double mass_ratio, double stell_luminosity_ratio,
               int mass_type, planet_pointer planet_list, double radius)
{
    double a, e, mass, crit_mass, eff_inner_bound, eff_outer_bound,
      planet_inner_bound, planet_outer_bound,
	  dust_inner_bound, dust_outer_bound,
      bound1, bound2, temp1;
	dust_pointer band;

	/*
	 *	Figure out the inner and outer limits at which a body can exist
	 *	about this body ('planet_inner_bound' and 'planet_outer_bound'):
	 */
    if (mass_type == MOON) {
        planet_head = NULL;
        planet_inner_bound = roche_limit(radius * 2.0);
    }
    else {
        planet_head = planet_list;
        planet_inner_bound = nearest_body(mass_ratio);
    }
    planet_outer_bound = farthest_body(mass_ratio);

	/*
	 *	Figure out the innermost and outermost extent of the dust/gas
	 *	cloud about the object.  The dust can't be any closer to the
	 *	primary than can be affected by a protoplanet with zero orbital
	 *	eccentricity at the minimum distance from the primary:
	 */
	if (mass_type == PLANET) {
		dust_inner_bound = inner_effect_limit(planet_inner_bound,
			0.0, PROTOPLANET_MASS);
		dust_outer_bound = stell_dust_limit(mass_ratio, 0.0, STAR);
		temp1 = outer_effect_limit(planet_outer_bound, 0.0, PROTOPLANET_MASS);
		if (dust_outer_bound > temp1) {
			dust_outer_bound = temp1;
		}
	}
	else if (mass_type == MOON) {
		dust_inner_bound = inner_effect_limit(planet_inner_bound,
			0.0, PROTOPLANET_MASS);
		dust_outer_bound = stell_dust_limit(mass_ratio, radius, PLANET);
		temp1 = outer_effect_limit(planet_outer_bound, 0.0, PROTOPLANET_MASS);
		if (dust_outer_bound > temp1) {
			dust_outer_bound = temp1;
		}
	}
	else {
		printf("ERROR: bad mass type in 'dist_masses'\n");
		exit(1);
	}
	/*
	 *  Set up a clean dust/gas cloud in a single band about the object:
	 */
    set_initial_conditions(dust_inner_bound, dust_outer_bound);

	/*
	 *	Inject proto-masses until all the dust about the central body
	 *	has been accumulated:
	 */
    while (dust_left) {
        e = random_eccentricity( );
        mass = PROTOPLANET_MASS;
#ifdef NOTUSED
        innermost_limit = inner_effect_limit(planet_inner_bound,
                             e,
                             mass);
        outermost_limit = outer_effect_limit(planet_outer_bound,
                             e,
                             mass);
        if (innermost_limit > inner_dust) {
            bound1 = planet_inner_bound;
        }
        else {
            bound1 = inner_dust;
        }
        if (outermost_limit < outer_dust) {
            bound2 = planet_outer_bound;
        }
        else {
            bound2 = outer_dust;
        }
        if ((bound2 - bound1) < 0.0) {
            return(NULL);
        }
#endif
		/*
		 *	Find the first dust/gas band with dust still present:
		 */
		band = dust_head;
		while ((band != NULL) && (band->dust_present == FALSE)) {
			band = band->next_band;
		}
		if (band == NULL) {
			printf("ERROR: dust band checking internal error\n");
			exit(1);
		}
		/*
		 *	Choose a location for the proto-mass that is somewhere
		 *	within gravitational effect range of the first band.  As
		 *	this is done for each band, the innermost band with dust
		 *	still remaining will move further and further from the primary
		 *	until all dust in the system has been accreted.
		 */
#ifdef NOTDEF
        bound1 = inner_effect_limit(band->inner_edge, e, mass);
        bound2 = outer_effect_limit(band->outer_edge, e, mass);
#endif
        bound1 = band->inner_edge;
        bound2 = band->outer_edge;
		if (bound1 < planet_inner_bound) {
			bound1 = planet_inner_bound;
		}
		if (bound2 > planet_outer_bound) {
			bound2 = planet_outer_bound;
		}
		if (planet_inner_bound > planet_outer_bound) {
			printf("ERROR: orbit bounding internal error\n");
			exit(1);
		}
        a = random_number(bound1, bound2);
        eff_inner_bound = inner_effect_limit(a, e, mass);
        eff_outer_bound = outer_effect_limit(a, e, mass);
        if (dust_available(eff_inner_bound, eff_outer_bound)) {
            if (flag_verbose >= LEVEL1) {
                if (mass_type == PLANET)
                    printf("  Injecting proto-planet (%4.2lg AU)\n", a);
                else printf("  Injecting proto-moon (%4.2lg AU)\n", a);
            }
            dust_density = DUST_DENSITY_COEFF * sqrt(mass_ratio)
                * exp(-ALPHA * power(a,(1.0 / N)));
			/*
			 *	Assume that dust is ten times more dense around
			 *	planets:
			 */
			if (mass_type == MOON) {
				dust_density = dust_density * 10.0;
			}
            crit_mass =critical_limit(a,e,stell_luminosity_ratio);
            mass = accrete_dust(mass,a,e,crit_mass);
            if ((mass != 0.0) && (mass != PROTOPLANET_MASS)) {
                coalesce_planetesimals(a,e,mass,crit_mass,
                               stell_luminosity_ratio,
                               mass_type);
            }
            else if (flag_verbose >= LEVEL2) {
                printf("    Neighbor too near (%lg AU).\n",a);
            }
        }
        else if (flag_verbose >= LEVEL2) {
            printf("    Not enough dust at %lg AU.\n",a);
        }
    }
    return(planet_head);
}

/*--------------------------------------------------------------------------*/
/*  This function checks if each planet is within the radius of the star    */
/*  or if it is at least close enough to be vaporized.  If either of these  */
/*  is the case, the planet will be deleted.                                */
/*  The head is the head of the list of planets, luminosity is a unitless   */
/*  ratio of the star's luminosity to that of the Sun's, and the radius of  */
/*  the star is given in AU.                                                */
/*--------------------------------------------------------------------------*/
planet_pointer check_planets(planet_pointer head, double luminosity, double star_radius)
{
    planet_pointer planet, absorbed_planet, vaporized_planet;
    double r_ecosphere, temperature;

    if (head == NULL)
        return(NULL);
    r_ecosphere = sqrt(luminosity);
    /*
     *  Start with the second planet on the list - the first one is
     *  always the primary star:
     */
    planet = head->next_planet;
    while (planet != NULL) {
        temperature = eff_temp(r_ecosphere, planet->a, ROCKY_AIRLESS_ALBEDO);
        if (planet->a <= star_radius) {
            /*
             *  The planet is inside the primary!  Zap it:
             */
            absorbed_planet = planet;
            planet = planet->next_planet;
            if (absorbed_planet == head)
                head = planet;
            free(absorbed_planet);
            if (flag_verbose >= LEVEL1)
                printf("  Planet absorbed by primary!\n");
        }
        else if (temperature >= 2000.0) {
            /*
             *  Too hot!  Zap it:
             */
            vaporized_planet = planet;
            planet = planet->next_planet;
            if (vaporized_planet == head)
                head = planet;
            free(vaporized_planet);
            if (flag_verbose >= LEVEL1)
                printf("  Planet vaporized by primary!\n");
        }
        else planet = planet->next_planet;
    }
    return(head);
}

/*--------------------------------------------------------------------------*/
/*  As input, this function receives a pointer to the head of the list of   */
/*  star structs.  It creates a corresponding planet struct for each of     */
/*  the star structs and returns a pointer to the head of the new planet    */
/*  list.                                                                   */
/*--------------------------------------------------------------------------*/
planet_pointer init_planet_list (star_pointer star_head)
{
    star_pointer star;
    planet_pointer planet;
    planet_pointer planet_list_head = NULL;

    for (star = star_head; star != NULL; star = star->next_star) {
        if ((planet = (planets *)malloc((unsigned)sizeof(planets))) == NULL) {
            perror("malloc'ing a new star/planet");
            exit(1);
        }
        if (flag_verbose >= LEVEL3) {
            printf("      Creating a new planet node for a star.\n");
        }
        planet->mass_type = STAR;
        planet->a = star->orbit_radius;
        planet->mass = star->stell_mass_ratio;
        planet->e = random_eccentricity();
        /*
         *  Insert the new planet in the planet list.  Keep the planet
         *  list sorted by distance from the primary:
         */
        planet_list_head = sorted_list_insert(planet_list_head, planet);
        /*
         *  Now create a double link between the corresponding star and
         *  planet node:
         */
        planet->star_ptr = star;
        star->planet_ptr = planet;
    }
    return(planet_list_head);
}

