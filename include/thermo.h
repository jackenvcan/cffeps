#ifndef THERMO_H
#define THERMO_H
/* NOTE ALL PARAMETERS PASSED TO AND RETURNED FROM 
   THESE THERMODYNAMIC ROUINES USE

      TEMPERATURE IN CELSIUS
      PRESSURE IN MILLIBARS
      HUMIDITY IN PERCENT  */

/* some thermodynamic constants */
static double 
       c_p = 1005.,
       g = 9.8,
       l_v = 2.5008E6,
       R_d = 287.05,
       chi_d = 0.286,
       epsilon = 0.622;

/* latent heat of vaporization [J/kg/K]
   This equation was derived as a linear interpolation of the data points in 
   Irabarne & Godson */
double l_v_Calc(double T /* oC */);

/* latent heat of fusion [J/kg/K]
   This equation was derived as a linear interpolation of the data points in 
   Irabarne & Godson */
double l_f_Calc(double T /* oC */);

/* latent heat of sublimation [J/kg/K] 
   This equation was supplied bydjboccip@cirrus.mit.edu (Dennis J. Boccippio)
   and was derived from the data points in Irabarne & Godson */
double l_s_Calc(double T /* oC */);

/* vapour pressure over water [mb] */
double e_w_Calc(double T /* oC */);

/* vapour pressure over ice [mb] */
double e_i_Calc(double T /* oC */);

/* The following equations are taken from Irabarne and Godson (1973)
   References are by chapter-equation */

/* mixing ratio [] */
double r_Calc(double p    /* mb */,
              double T_d  /* oC */);

/* dew point temperature [oC] */
double T_d_Calc(double T  /* oC */, 
                double rh /* % */);

/* wet-bulb temperature [oC] */
double T_w_Calc(double p   /* mb */, 
                double T   /* oC */, 
                double T_d /* oC */);

/* dry adiabat line */
double adiabat(double T,    /* oC */
               double theta /* K */);

/* Variation of T with height */
double T_wrt_z(double T,  /* oC */
               double from_z, /* m */
               double to_z  /* m */);

/* Variation of rh with height */
double rh_wrt_z(double T,  /* oC */
                double rh, /* % */
                double from_z, /* m */
                double to_z  /* m */);
#endif
