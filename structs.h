/*----------------------------------------------------------------------*/
/*                             structs.h                                */
/*                                                                      */
/*  Header file containing all structures used in starform.             */
/*----------------------------------------------------------------------*/
/* (c) Copyright Matt Burdick 1991 - All Rights Reserved                */
/*----------------------------------------------------------------------*/
typedef struct dust_struct  *dust_pointer;
typedef struct planets_struct  *planet_pointer;
typedef struct star_struct *star_pointer;
typedef struct sys_struct *sys_pointer;

typedef struct sys_struct {
	star_pointer primary_star;
	planet_pointer inner_planet;
} star_system;

typedef struct star_struct {
	double orbit_radius;                       /* in AU */
	double age;                                /* in years */
	double stell_mass_ratio;                   /* in solar masses */
	char lum_id;                               /* command-line luminosity */
	char spec_class;                           /* command-line spec class */
	int spec_num;                              /* command-line spec number */
	char star_type[CLASSIFICATION_SIZE];  /* stellar classification	    */
	int lum_type;
	double stell_luminosity_ratio;
	double stell_radius;                       /* in AU */
	double main_seq_life;                      /* in years */
	double r_ecosphere;                        /* in AU */
	double r_greenhouse;
	star_pointer next_star;
	planet_pointer planet_ptr;
} stars;

typedef struct planets_struct {
	double a;		/* semi-major axis of the orbit (in AU)*/
	double e;		/* eccentricity of the orbit	     */
	double mass;		/* mass (in solar masses)	     */
	int mass_type;		/* indicates star, planet, moon, etc */
	int orbit_zone;         /* the 'zone' of the planet          */
	double radius;		/* equatorial radius (in km)	     */
	double density;		/* density (in g/cc)		     */
	double orb_period;   	/* length of the local year (days)   */
	double day;		/* length of the local day (hours)   */
	int resonant_period;	/* TRUE if in resonant rotation   */
	int axial_tilt;		/* units of degrees		     */
	double esc_velocity;	/* units of cm/sec		     */
	double surf_accel;  	/* units of cm/sec2		     */
	double surf_grav;   	/* units of Earth gravities	     */
	double rms_velocity;	/* units of cm/sec		     */
	double molec_weight;	/* smallest molecular weight retained*/
	double volatile_gas_inventory;
	double surf_pressure;	/* units of millibars (mb)	     */
	int greenhouse_effect;	/* runaway greenhouse effect?	*/
	double boil_point;	/* the boiling point of water (Kelvin)*/
	double albedo;		/* albedo of the planet		     */
	double surf_temp;   	/* surface temperature in Kelvin     */
	double hydrosphere;	/* fraction of surface covered	     */
	double cloud_cover;	/* fraction of surface covered	     */
	double ice_cover;	/* fraction of surface covered	     */
	planet_pointer first_moon;
	planet_pointer next_planet;
	star_pointer star_ptr;
} planets;


typedef struct dust_struct {
	double inner_edge;
	double outer_edge;
	int dust_present;
	int gas_present;
	dust_pointer next_band;
} dust;

typedef struct Spectral_Info {
	char spec_class;
	int spec_num;
	double max_mass;
	int percentage;
} spectral_info;

