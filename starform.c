/*----------------------------------------------------------------------*/
/*                            starform.c                                */
/*                                                                      */
/* Main routine for running the star system accretion process and then  */
/* visiting each planet to determine characteristics.                   */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
/*
 *    C language includes:
 */
#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#include    <sys/types.h>
#include    <sys/timeb.h>

#ifdef MSDOS
#include	<string.h>
#include    <stddef.h>
#include    <malloc.h>
#include    <stdlib.h>
#include    <float.h>
#endif

#ifdef BSD
#include	<strings.h>
#else
#include	<string.h>
#endif

/*
 *    Starform includes:
 */
#include    "config.h"
#include    "const.h"
#include    "structs.h"

#ifdef MSDOS
#include    "protos.h"
#endif

/*
 *    Global variables used during planetary accretion:
 */
star_system sys;
int resonance;


/*
 *    Global variables used for command-line parameters:
 */
unsigned flag_seed =    FALSE;
unsigned flag_verbose = FALSE;
int flag_graphic =      FALSE;
int flag_moons =        FALSE;
int flag_startype =     FALSE;
int flag_tec =          FALSE;  /* for Dave Allen's "tec" program */

/*
 *    Externals from C library not elsewhere declared:
 */
extern void srand();

/*
 *    Externals from other starform files (mainly enviro.c):
 */
extern double random_number();
extern double luminosity();
extern double star_radius();
extern char *classify();
extern double star_mass();
extern double star_age();
extern planet_pointer dist_masses();
extern planet_pointer check_planets();
extern planet_pointer init_planet_list();
extern int orb_zone();
extern double empirical_density();
extern double volume_density();
extern double volume_radius();
extern double kothari_radius();
extern double period();
extern double day_length();
extern int inclination();
extern double escape_vel();
extern double accel();
extern double rms_vel();
extern double molecule_limit();
extern double about();
extern int grnhouse();
extern double gravity();
extern double vol_inventory();
extern double pressure();
extern double boiling_point();
extern int  verify_startype();
extern void startype_error();
extern int rand_type();
extern double rand_star_mass();
extern void display_system();
extern void iterate_surface_temp();

void usage();
void init();
void generate_stellar_system();


/*--------------------------------------------------------------------------*/
/*   The main function decodes all of the command-line parameters, then     */
/*   calls 'generate_stellar_system' and 'display_system'.  Currently, the  */
/*   -g flag doesn't work.  Eventually, it will provide some sort of        */
/*   graphical output, but I can't decide whether MS-Windows or X windows   */
/*   would be better (after all, I don't have an X windows workstation at   */
/*   home).                                                                 */
/*--------------------------------------------------------------------------*/
void main (argc, argv)
int argc;
char *argv[];
{
    char *c, *progname;
    int skip;
    int errornum;
    int radius;
    star_pointer *star = &sys.primary_star;

/*
 *    Grab all the command-line parameters:
 */
    progname = argv[0];
    while (--argc > 0) {
        if ((*++argv)[0] != '-')
            usage(progname);
        for (c = argv[0]+1, skip=FALSE;
             (*c != '\0') && (!(skip));
             c++)
            switch (*c) {
            case 'g':    /* display graphically */
                ++flag_graphic;
                break;
            case 'm':    /* generate moons for planets */
                ++flag_moons;
                break;
            case 's':    /* set random seed */
                flag_seed = (unsigned) atoi(&(*++c));
                skip = TRUE;
                break;
            case 'v':    /* increment verbosity */
                flag_verbose = (unsigned) atoi(&(*++c));
                skip = TRUE;
                break;
            case 'T':    /* generate "tec" output */
                ++flag_tec;
                break;
            case 't': /* specify star type */
                ++flag_startype;
                /*
                 *  Allocate memory for the new star in the linked list:
                 */
                if ((*star = (stars *)malloc((unsigned)sizeof(stars))) == NULL) {
                    perror("malloc'ing memory for a star");
                    exit(1);
                }
                (*star)->next_star = NULL;
                /*
                 *  Figure out stellar class and luminosity information:
                 */
                if (sscanf(++c, "%c%d%c/%d", &(*star)->spec_class,
                                      &(*star)->spec_num,
                                      &(*star)->lum_id,
                                      &radius) != 4) {
                    usage(progname);
                }
                (*star)->orbit_radius = (double)radius;
                if ((errornum = verify_startype((*star)->lum_id,
                                                (*star)->spec_num,
                                                (*star)->spec_class)) != 0) {
                    startype_error(errornum, (*star)->spec_class,
                                              (*star)->spec_num,
                                             (*star)->lum_id);
                    usage(progname);
                }
                /*
                 *  Translate the command-line 'lum_id' into a 'lum_type':
                 */
                switch ((*star)->lum_id) {
                    case 'S':
                        (*star)->lum_type = SUPERGIANT;
                        break;
                    case 'G':
                        (*star)->lum_type = GIANT;
                        break;
                    case 'D':
                        (*star)->lum_type = WHITE_DWARF;
                        break;
                    case 'M':
                    default:
                        (*star)->lum_type = MAIN_SEQUENCE;
                        break;
                }
                skip = TRUE;
                star = &((*star)->next_star);
                break;
            default:
            case '?':
                usage(progname);
            }
    }
/*
 *    Now do all the hard work:
 */
    init();
    generate_stellar_system();
    display_system(&sys);
}

/*--------------------------------------------------------------------------*/
/*   Tell the user what kind of command-line parameters are possible.       */
/*--------------------------------------------------------------------------*/
void usage (progname)
char *progname;
{

    fprintf(stderr,
        "%s: Usage: %s [-g] [-m] [-s#] [-v#] [-tl#l/#]\n", progname);
    fprintf(stderr,
        "\t -g        Display graphically (unimplemented)\n");
    fprintf(stderr,
        "\t -m        Generate moons for each planet\n");
    fprintf(stderr,
        "\t -s#       Use # as the seed for random number generation\n");
    fprintf(stderr,
        "\t -v#       Set the verbosity level to # (default is 0)\n");
    fprintf(stderr,
        "\t -tl#l/#   Choose the spectral type, luminosity class, and orbit\n");
    exit (1);
}

/*--------------------------------------------------------------------------*/
/*   Initialize the random-number generator.                                */
/*--------------------------------------------------------------------------*/
void init()
{
    unsigned seed;
    struct timeb grap;

    if (flag_seed)
        seed = flag_seed;
    else {
        ftime(&grap);
        seed = (unsigned)((grap.time%100000)+grap.millitm);
    }
    (void)srand(seed);
    printf("Starform - V%s\n", VERSION);
    printf("Random number seed - %u\n", seed);
}

/*--------------------------------------------------------------------------*/
/*   First, find out what kind of stars are in this system, then use the    */
/*   'dist_masses' function to accrete dust and gasses into planets.        */
/*   Finally, loop through each planet finding the physical                 */
/*   characteristics of each one.                                           */
/*--------------------------------------------------------------------------*/
void generate_stellar_system()
{
    planet_pointer planet;
    planet_pointer moon;
    star_pointer star;
    star_pointer previous_star;
    char *buf;
    int temp;    /* Used in calculating the number of stars in a system */
    int star_number;/* The number of stars in this system.                 */

    /*
     *  Build up the list of stars in this system.  If the user specified
     *  the stars on the command-line, use those.  Otherwise, randomly
     *  determine how many stars to generate, then create star types
     *  and orbital distances for each of those.
     */
    if (flag_startype) {
        for (star = sys.primary_star; star != NULL; star = star->next_star) {
            star->stell_mass_ratio = star_mass(star->lum_type,
                                               star->spec_class,
                                               star->spec_num);
            if (star->stell_mass_ratio == 0.0) {
                fprintf(stderr,"ERROR: white dwarfs are rarely type M\n");
                exit(1);
            }
            /*
             *  Create the normal text description of the spectral class:
             */
            switch (star->lum_type) {
                case GIANT:
                    sprintf(star->star_type, "%c%d III\0", star->spec_class, star->spec_num);
                    break;
                case SUPERGIANT:
                    sprintf(star->star_type, "%c%d Ia\0", star->spec_class, star->spec_num);
                    break;
                case WHITE_DWARF:
                    sprintf(star->star_type, "D%c%d\0", star->spec_class, star->spec_num);
                    break;
                case MAIN_SEQUENCE:
                default:
                    sprintf(star->star_type, "%c%d V\0", star->spec_class, star->spec_num);
                    break;
            }
        }
    }
    else {
        /*
         *  Decide how many stars should be in this system.  The percentage
         *  of double, triple, and quadruple star systems is basically pulled
         *  from a hat - the best estimates of the actual frequencies of
         *  these kind of systems I could find said only that "more than
         *  half of all stars are members of multiple star systems".
         */
        temp = (int)random_number(1.0, 100.0);
        if (temp <= 45) {
            star_number = 1;
        }
        else if (temp <= 80) {
            star_number = 2;
        }
        else if (temp <= 95) {
            star_number = 3;
        }
        else {
            star_number = 4;
        }
        if (flag_verbose >= LEVEL1) {
            printf("  Creating system with %d stars.\n", star_number);
        }
        /*
         *  Determine basic characteristics of all the stars:
         */
        for (temp = 1; temp <= star_number; temp++) {
            if ((star = (stars *)malloc((unsigned)sizeof(stars))) == NULL) {
                perror("malloc'ing memory for a star");
                exit(1);
            }
            if (temp == 1) {    /* Is this the first star? */
                sys.primary_star = star;
                star->orbit_radius = 0.0;
            }
            else {
                previous_star->next_star = star;
                star->orbit_radius = random_number(1.0, 150.0);
            }
            previous_star = star;
            star->lum_type = rand_type();
            star->stell_mass_ratio = rand_star_mass(star->lum_type);
            buf = classify(star->stell_mass_ratio, star->lum_type);
            (void)strncpy(star->star_type, buf, CLASSIFICATION_SIZE);
            star->next_star = NULL;
        }
    }
    /*
     *  The rest of the stellar characteristics depend on those above that
     *  have either been specified on the command line or generated
     *  randomly.
     */
    for (star = sys.primary_star; star != NULL; star = star->next_star) {
        star->stell_luminosity_ratio = luminosity(star->stell_mass_ratio,
        star->lum_type);
        if ((star->star_type[0] == 'K') || (star->star_type[0] == 'M')) {
            star->stell_radius = star_radius(star->stell_mass_ratio, star->lum_type,
            TRUE);
        }
        else {
            star->stell_radius = star_radius(star->stell_mass_ratio, star->lum_type,
            FALSE);
        }
        star->main_seq_life = 1.1E10 * (star->stell_mass_ratio
            / star->stell_luminosity_ratio);
        if (star->main_seq_life < 1.0E6) {
            star->main_seq_life = 1.0E6;
        }
        star->age = star_age(star->main_seq_life);
        star->r_ecosphere = sqrt(star->stell_luminosity_ratio);
        star->r_greenhouse = star->r_ecosphere * GREENHOUSE_EFFECT_CONST;
    }
    if (flag_verbose >= LEVEL1) {
        printf("  Begin building main planetary orbits:\n");
    }
/*
 *  Now that we have the star information, build a planetary system
 *  through accretion.  Start by adding all the stars into the planet
 *  list, then use 'dist_masses' to inject protoplanets until there's
 *  no more gas or dust to collect:
 */
    sys.inner_planet = init_planet_list(sys.primary_star);
    sys.inner_planet = dist_masses(sys.primary_star->stell_mass_ratio,
                       sys.primary_star->stell_luminosity_ratio,
                       PLANET, sys.inner_planet, 0.0);
/*
 *  Now check if each planet is within the radius of the primary star or
 *  at least close enough to be vaporized:
 */
    sys.inner_planet = check_planets(sys.inner_planet,
            sys.primary_star->stell_luminosity_ratio,
            sys.primary_star->stell_radius);
    if (flag_verbose >= LEVEL1) {
        printf("  Finished building planetary orbits\n");
    }
    for (planet=sys.inner_planet;
         planet != NULL;
         planet = planet->next_planet) {
        /*
         *  If this 'planet' is really a star, skip it:
         */
        if (planet->mass_type == STAR) {
            continue;
        }
        planet->orbit_zone =orb_zone(planet->a,
                         sys.primary_star->stell_luminosity_ratio);
        if (planet->mass_type == GAS_GIANT) {
            planet->density = empirical_density(planet->mass,
                                planet->a,
                                planet->mass_type,
                                sys.primary_star->stell_luminosity_ratio);
            planet->radius = volume_radius(planet->mass,
                               planet->density);
        }
        else {
            planet->radius = kothari_radius(planet->mass,
                            planet->mass_type,
                            planet->orbit_zone);
            planet->density = volume_density(planet->mass,
                             planet->radius);
        }
/*
 *  Build the planet's moons if moons were specified on the command line
 *  and the 'planet' isn't really a companion star:
 */
        if (flag_moons && (planet->mass_type != STAR)) {
            planet->first_moon =
                dist_masses(planet->mass,
                        sys.primary_star->stell_luminosity_ratio,
                        MOON,
                        NULL,
                        planet->radius);
            if (flag_verbose >= LEVEL1) {
                printf("  Built moon orbits for a planet\n");
            }
            for (moon=planet->first_moon;
                moon != NULL;
                moon = moon->next_planet) {
                if (moon->mass_type == GAS_GIANT) {
                    moon->density = empirical_density(moon->mass,
                                        planet->a,
                                        moon->mass_type,
                                        sys.primary_star->r_ecosphere);
                    moon->radius = volume_radius(moon->mass,
                                      moon->density);
                }
                else {
                    moon->radius = kothari_radius(moon->mass,
                                    moon->mass_type,
                                    planet->orbit_zone);
                    moon->density = volume_density(moon->mass,
                                     moon->radius);
                }
                moon->surf_accel = accel(moon->mass, moon->radius);
                moon->surf_grav = gravity(moon->surf_accel);
            }
        }
        else {
            planet->first_moon = NULL;
        }
        planet->orb_period = period(planet->a,
                        planet->mass,
                        sys.primary_star->stell_mass_ratio);
        planet->day = day_length(planet->mass,
                     planet->radius,
                     planet->e,
                     planet->density,
                     planet->a,
                     planet->orb_period,
                     planet->mass_type,
                     sys.primary_star->stell_mass_ratio,
                     sys.primary_star->age);
        planet->resonant_period = resonance;
        planet->axial_tilt = inclination(planet->a);
        planet->esc_velocity = escape_vel(planet->mass,
                          planet->radius);
        planet->surf_accel = accel(planet->mass,planet->radius);
        planet->rms_velocity = rms_vel(MOL_NITROGEN,planet->a,
             sys.primary_star->stell_luminosity_ratio);
        planet->molec_weight = molecule_limit(planet->mass,
                              planet->radius);
        if (planet->mass_type == GAS_GIANT) {
            planet->surf_grav = 0.0;
            planet->greenhouse_effect = FALSE;
            planet->volatile_gas_inventory = 0.0;
            planet->surf_pressure = 0.0;
            planet->boil_point = 0.0;
            planet->hydrosphere = 0.0;
            planet->albedo = about(GAS_GIANT_ALBEDO,0.1);
            planet->surf_temp = 0.0;
        }
        else {
            planet->surf_grav = gravity(planet->surf_accel);
            planet->greenhouse_effect=grnhouse(planet->orbit_zone,
                               planet->a,
                               sys.primary_star->r_greenhouse);
            planet->volatile_gas_inventory =
                vol_inventory(planet->mass,
                          planet->esc_velocity,
                          planet->rms_velocity,
                          sys.primary_star->stell_mass_ratio,
                          planet->orbit_zone,
                          planet->greenhouse_effect);
            planet->surf_pressure = pressure(planet->volatile_gas_inventory,
                     planet->radius,
                     planet->surf_grav);
            if (planet->surf_pressure == 0.0)
                planet->boil_point = 0.0;
            else planet->boil_point = boiling_point(planet->surf_pressure);
            iterate_surface_temp(&(planet), sys.primary_star->r_ecosphere);
        }
    }
}


