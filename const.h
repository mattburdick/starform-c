/*----------------------------------------------------------------------*/
/*                               const.h                                */
/*                                                                      */
/*  Header file containing all the nasty constants used in starform.    */
/*  For fun, try tweaking ECCENTRICITY_COEFF (determining the           */
/*  eccentricity of your planet's orbits about the star) and            */
/*  DUST_DENSITY_COEFF (determining the density of the initial dust     */
/*  cloud about the protostar) to see systems with different mass       */
/*  distribution characteristics.                                       */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
#ifndef RAND_MAX
#ifdef LONG_RAND
#define RAND_MAX        (2147483647.0)
#else
#define RAND_MAX        (32767.0)
#endif /* ifdef LONG_RAND */
#endif /* ifndef RAND_MAX */

#define BAD_LUMINOSITY	1
#define BAD_SPECTRA		2
#define BAD_MOD		4

#define PI                      (3.1415926536)
#define	RADIANS_PER_ROTATION	(2.0 * PI)
#define TRUE                    (1)
#define FALSE                   (0)
#define ECCENTRICITY_COEFF      (0.077)		/* Dole's was 0.077         */
#define CLOUD_ECCENTRICITY      (0.2)
#define	PROTOPLANET_MASS	(1.0E-15)	/* Units of solar masses    */
#define CHANGE_IN_EARTH_ANG_VEL (-1.3E-15)	/* Units of radians/sec/year*/
#define SOLAR_MASS_IN_GRAMS	(1.989E33)	/* Units of grams           */
#define EARTH_MASS_IN_GRAMS	(5.977E27)	/* Units of grams           */
#define MOON_MASS_IN_GRAMS	(7.47E25)	/* Units of grams           */
#define EARTH_RADIUS		(6.378E8)	/* Units of cm		    */
#define EARTH_DENSITY		(5.52)		/* Units of g/cc	    */
#define KM_EARTH_RADIUS		(6378.0)	/* Units of km              */
#define EARTH_ACCELERATION	(981.0)		/* Units of cm/sec2         */
#define EARTH_AXIAL_TILT	(23.4)		/* Units of degrees         */
#define EARTH_EXOSPHERE_TEMP	(1273.0)	/* Units of degrees Kelvin  */
#define SUN_MASS_IN_EARTH_MASSES (332775.64)
#define SOLAR_RADII_PER_AU    (4.652E-3) /* Number of solar radii per AU */
#define EARTH_EFFECTIVE_TEMP	(255.0)		/* Units of degrees Kelvin  */
#define EARTH_ALBEDO		(0.3)
#define	CLOUD_COVERAGE_FACTOR	(1.839E-8)	/* Km2/kg                   */
#define EARTH_WATER_MASS_PER_AREA	(3.83E15)/* grams per square km     */
#define	EARTH_SURF_PRES_IN_MILLIBARS	 (1000.0)
#define EARTH_CONVECTION_FACTOR	(0.43)		/* from Hart, eq.20         */
#define FREEZING_POINT_OF_WATER	(273.0)		/* Units of degrees Kelvin  */
#define DAYS_IN_A_YEAR		(365.256)	/* Earth days per Earth year*/
#define GAS_RETENTION_THRESHOLD	(6.0)		/* ratio of esc vel to RMS vel */
#define GAS_GIANT_ALBEDO	(0.5)		/* albedo of a gas giant    */
#define CLOUD_ALBEDO		(0.52)
#define ROCKY_AIRLESS_ALBEDO	(0.07)
#define ROCKY_ALBEDO		(0.15)
#define WATER_ALBEDO		(0.04)
#define AIRLESS_ICE_ALBEDO	(0.5)
#define ICE_ALBEDO		(0.7)
#define SECONDS_PER_HOUR	(3600.0)
#define CM_PER_AU		(1.495978707E13)/* number of cm in an AU    */
#define CM_PER_KM		(1.0E5)		/* number of cm in a km     */
#define KM_PER_AU               (CM_PER_AU / CM_PER_KM)
#define CM_PER_METER		(100.0)
#define MILLIBARS_PER_BAR	(1000.0)
#define KELVIN_CELCIUS_DIFFERENCE	(273.0)
#define GRAV_CONSTANT		(6.672E-8)	/* units of dyne cm2/gram2  */
#define GREENHOUSE_EFFECT_CONST	(0.93)		/* affects inner radius..   */
#define MOLAR_GAS_CONST		(8314.41)	/* units: g*m2/(sec2*K*mol) */
#define K			(50.0)		/* K = gas/dust ratio       */
#define B			(1.2E-5)	/* Used in Crit_mass calc   */
#define DUST_DENSITY_COEFF	(6.0E-3)	/* A in Dole's paper        */
#define	ALPHA			(5.0)		/* Used in density calcs    */
#define	N			(3.0)		/* Used in density calcs    */
#define J			(1.46E-19)	/* Used in day-length calcs (cm2/sec2 g) */
#define TRIVIAL_MASS            (1.0E-14)       /* Units of solar masses    */
#define TEMP_ITERATION_LIMIT	(101)		/* Limit on temp iterations */
#define CLASSIFICATION_SIZE	(10)		/* Size of star_type field  */


/*  Now for a few molecular weights (used for RMS velocity calcs):     */
/*  This table is from Dole's book "Habitable Planets for Man", p. 38  */

#define ATOMIC_HYDROGEN		(1.0)	/* H   */
#define MOL_HYDROGEN		(2.0)	/* H2  */
#define HELIUM			(4.0)	/* He  */
#define ATOMIC_NITROGEN		(14.0)	/* N   */
#define ATOMIC_OXYGEN		(16.0)	/* O   */
#define METHANE			(16.0)	/* CH4 */
#define AMMONIA			(17.0)	/* NH3 */
#define WATER_VAPOR		(18.0)	/* H2O */
#define NEON			(20.2)	/* Ne  */
#define MOL_NITROGEN		(28.0)	/* N2  */
#define CARBON_MONOXIDE		(28.0)	/* CO  */
#define NITRIC_OXIDE		(30.0)	/* NO  */
#define MOL_OXYGEN		(32.0)	/* O2  */
#define HYDROGEN_SULPHIDE	(34.1)	/* H2S */
#define ARGON			(39.9)	/* Ar  */
#define CARBON_DIOXIDE		(44.0)	/* CO2 */
#define NITROUS_OXIDE		(44.0)	/* N2O */
#define NITROGEN_DIOXIDE	(46.0)	/* NO2 */
#define OZONE			(48.0)	/* O3  */
#define SULPH_DIOXIDE		(64.1)	/* SO2 */
#define SULPH_TRIOXIDE		(80.1)	/* SO3 */
#define KRYPTON			(83.8)	/* Kr  */
#define XENON			(131.3)	/* Xe  */

/*
 *  Define's for indicating body type:
 */
#define STAR       0
#define PLANET     1
#define GAS_GIANT  2
#define MOON       3

/*
 *  The following are for use in star generation in stars.c:
 */
#define GIANT			1
#define SUPERGIANT		2
#define WHITE_DWARF		3
#define MAIN_SEQUENCE 	4

/*
 *  The following defines are used in the kothari_radius function in
 *  file enviro.c.
 */
#define	A1_20			(6.485E12)	/* All units are in cgs system.  */
#define A2_20			(4.0032E-8)	/*   ie: cm, g, dynes, etc.      */
#define BETA_20			(5.71E12)


/*
 *  The following defines are used in determining the fraction of a planet
 *  covered with clouds in function cloud_fraction in file enviro.c.
 */
#define	Q1_36			(1.258E19)	/* grams    */
#define Q2_36			(0.0698)	/* 1/Kelvin */


/*
 *  The following defines are for verbosity:
 */
#define  LEVEL1  1
#define  LEVEL2  2
#define  LEVEL3  3
#define  LEVEL4  4


