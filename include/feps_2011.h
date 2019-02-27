#ifndef FEPS_2011_H
#define FEPS_2011_H

#define MAX_LEVELS 100
#define MAX_SPECIES 20
#define MAX_EMODELS 10
#define SPECIES_SIZE 25
#define MAX_HOURS 72
#define MAX_CWFIS 30
#define NB_CHAR_PATH 240
#define MAX_TIMESTEPS 1440 // 24 hours of 1 minute timesteps will be the maximum array size


/* FEPS structure 
   much in common with FBP */

typedef struct
{
/* inputs */
   char method[32],	   /* method of employment: AB= Alberta Plume Study; ICAO=Standard Atmosphere test; anything else=cmc calculations */
		profilename[NB_CHAR_PATH],/* name of the WRF->Hysplit profile filename (profile0_.txt) */
		csvfilename[NB_CHAR_PATH],/* name of Peter's hotspot file name */
		emissionsfilename[NB_CHAR_PATH],/* name of emissions factor matrix file name */
		emissionsmodelsfilename[NB_CHAR_PATH],/* name of emissions models file name */
        outfilename[NB_CHAR_PATH],/* name of the output file name */
        previous[NB_CHAR_PATH],   /* name of the previous time step's ini file (required to ascertain Qo) */
        shape[32],       /* shape of entrainment cloud: line/wedge, crescent, ellipse */
        type[32],		/* type of upper air profile method used (average lapse rate, dry, wet) */
        ID[64];          /* fire ID */
  
   int header;    /* print header info 1=Yes, 0=No */
   int ddate;     /* detection date of fire YYYYMMDD */
   int dtime;     /* detection time of fire HHMM (time all values are assumed to be collected) */
   int fdate;	  /* forecast (valid) date YYYYMMDD */
   int ftime;	  /* forecast (valid) time HHMM */
   int timezone;  /* time zone - hour offset from UTC [standard time] */
   int obs;       /* the time of the plume observation use in the Alberta study (HHMM) */
   int LDT;       /* local daylight time: 1=Yes, 0=No */
   int diurnal;   /* Diurnal adjustment: 1=Yes, 0=No */
   int sinks;     /* include the sinks terms: 1=Yes, 0=No */
   int radiation; /* radiation method: < 0 Byram's 1200/8600; =0 A_wall/A_top; > 0  value as % of Qfire */
   int print;     /* print: 1=Yes, 0=No */
   int newfire;   /* whether the hotpsot is a new fire (start fresh) or old (use existing energy and emissions data) 1=Yes, 0=No */
   double cmcUTCo;	  /* hour of first fire report (0...72) */

   double temp,        /* surface temperature (oC) */
	      rh,              /* humidity (%) */
          ws,              /* wind speed (km/h) */
	      precip, 		   /* precipitation (mm) */
          DMC,
          DC,
          pressure,        /* MSL pressure (mb) */
         
          growth,          /* change in fire size over last timestep (ha) */
          residencetime,   /* period used to calculate fire size of base of plume (hours) */
          residencegrowth, /* fire growth over residence time (hours) */
          area,            /* size of fire at dtime (ha) */
	      estarea,		   /* estimated fire size at time of detection (ha) */
	  
          maxarea,         /* maximum growth size (used in Alberta study) (ha) */
          perimeter,       /* perimeter at dtime (m) */

	      reset,		   /* reset time of smoke plume calculations [decimal hours LST] (off if < 0) */
	      thstart,		   /* time to start (exclusive) top-hat fire growth [decimal hours LST] */
	      thend,		   /* time to end (inclusive) top-hat fire growth [decimal hours LST] */
         
          lapse,           /* lapse rate (oC/m) */

          ts,              /* temperature at surface (oC) */
          t850,            /* temperature at 850 mb */
          t700,            /* temperature at 700 mb */
          t500,            /* temperature at 500 mb */
          t250,            /* temperature at 250 mb */
          zs,              /* height at surface (m) */
          z850,            /* height at 850 mb */
          z700,            /* height at 700 mb */
          z500,            /* height at 500 mb */
          z250,            /* height at 250 mb */
          zt,              /* height of the tropopause (essentially the maximum height */
          
          alpha,           /* entrainment half-angle (o) */
          timestep,        /* timestep for modelling black body radiation loss (decimal hours) */
          elapsed,		   /* elapsed time (decimal hours) */
          
          Qo,              /* amount of energy previously injected into atmosphere (J) */
          Qfire,           /* energy of the fire: used when directly calculating the plume height from energy (in stand alone module) */
          Qplume,          /* energy in the plume: used when directly calculating the plume height from energy (in stand alone module) */
          Qf2Qp,           /* conversion factor from fire to plume: used when directly calculating the plume height from energy (in stand alone module) */
		  
          Qflaming,		   /* Qplume broken down into flaming, smoldering and residual components (used in emissions.c) */
          Qsmoldering,     /* Qplume broken down into flaming, smoldering and residual components (used in emissions.c) */
          Qresidual,       /* Qplume broken down into flaming, smoldering and residual components (used in emissions.c) */
	  
	      Qs[MAX_TIMESTEPS+1],

          dz,              /* plume height (m) */
	      dzmin,		   /* minimum plume height (m)*/
          M,               /* mass of the plume volume *(tonnes) */

          Flaming,         /* consumption during the flaming stage [kg/m2] (used in emissions.c) */
          Smoldering,      /* consumption during the smoldering stage [kg/m2] (used in emissions.c) */
          Residual,        /* consumption during the residual stage [kg/m2] (used in emissions.c) */
	      totalemissions,  /* total emissions from the fire (tonnes) (used in emissions.c) */
	      r_smoke,		   /* mixing ratio of smoke to clear air (g/kg) */
      
          DOB,             /* depth of burn (used in emissions.c) */
          DOB_F,           /* depth of burn - flaming (used in emissions.c) */
          DOB_S;           /* depth of burn - smouldering (used in emissions.c) */
         
}  FEPS;

typedef struct                        /* emissions data*/
{
   char	*species[MAX_SPECIES],
            emissionsfilename[NB_CHAR_PATH];/* name of emissions factor matrix file name */

   double a[MAX_EMODELS][MAX_SPECIES],	/* emissions factor [g/kg] of flaming combustion */
          b[MAX_EMODELS][MAX_SPECIES],	/* emissions factor [g/kg] of smoldering combustion */
          c[MAX_EMODELS][MAX_SPECIES],	/* emissions factor [g/kg] of residual combustion */

         f[MAX_SPECIES][MAX_TIMESTEPS+1],   /* cumulative emissions for flaming combustion */
         s[MAX_SPECIES][MAX_TIMESTEPS+1],   /* cumulative emissions for smoldering combustion */
         r[MAX_SPECIES][MAX_TIMESTEPS+1],   /* cumulative emissions for residual combustion */
/* we need to consider both end points UTC(t=0) and UTC(t=MAX_TIMESTEPS) */

         residence_flaming,    /* residence time for flaming combustion (hours) */
         residence_smoldering, /* residence time for smoldering combustion (hours)  - follows flaming */
         residence_residual;   /* residence time for residual combustion (hours)  - follows smoldering */
         
      int nModels, nSpecies,     /* number of models and species */
         emodel[MAX_EMODELS];  /* this allows for different emissions models based on CWFIS fuel types */

} EMISSIONS;

typedef struct
{
   char	*FuelType[MAX_CWFIS];
   int     nModels,                        /* number of models  (-1 = generic model) */
           emodel[MAX_CWFIS];              /* this allows for different emissions models based on CWFIS fuel types */
        
} EMODELS;

typedef struct
{
	double P[MAX_LEVELS],
		   Z[MAX_LEVELS],
		   T[MAX_LEVELS],
		   Td[MAX_LEVELS],
		   W[MAX_LEVELS],
		   D[MAX_LEVELS];
    int iLevels;
} UA;

#endif
