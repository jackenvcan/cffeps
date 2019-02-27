#include <stdio.h>
#include <math.h>
#include "thermo.h"

/* NOTE ALL PARAMETERS PASSED TO AND RETURNED FROM 
   THESE THERMODYNAMIC ROUINES USE

      TEMPERATURE IN CELSIUS
      PRESSURE IN MILLIBARS
      HUMIDITY IN PERCENT  */


/*** latent heat of vaporization [J/kg/K] ***********************************/
/* This equation was derived as a linear interpolation of the data points in 
   Irabarne & Godson */

double l_v_Calc(double T /* oC */)
{
   if (T > -273.16)
      return(1E6*(3.176931 - 0.002465676 * (T+273.16)));
   else
      return (-999.);
}

/*** latent heat of fusion [J/kg/K] *****************************************/
/* This equation was derived as a linear interpolation of the data points in 
   Irabarne & Godson */

double l_f_Calc(double T /* oC */)
{
   if (T > -273.16)
      return(1E6*(-0.3685415 + 0.002584857 * (T+273.16)));
   else
      return (-999.);
}

/*** latent heat of sublimation [J/kg/K] *************************************/
/* This equation was supplied bydjboccip@cirrus.mit.edu (Dennis J. Boccippio)
   and was derived from the data points in Irabarne & Godson */

double l_s_Calc(double T /* oC */)
{
   if (T > -273.16)
      return(1E6*(2.637 + 0.0017 * (T+273.16) - 3.5629E-6 * (T+273.16)*(T+273.16)));
   else
      return (-999.);
}

/*** vapour pressure over water [mb] ****************************************/
double e_w_Calc(double T /* oC */)
{
double vp;

    if (T > -273.16)
    {

/* Alternative vapour pressure equations 
            vp = pow(10, 9.4041 - 2354/(T+273.16));                                    / * IV - 52 * /
            vp = pow(10, -2937.4/(T+273.16) - 4.9283*log10(T+273.16) + 23.5470);    / * IV - 53 * /
            vp = 6.1078*exp(17.269*T/(T+237.3));                                    / * Teten */


/* Goff-Gratch formula (List 1951) */
            vp = pow(10,-7.90298*(373.16/(T+273.16)-1.)
                        +5.02808*log10(373.16/(T+273.16))
                        -1.3816E-07*(pow(10,11.344*(1.-((T+273.16)/373.16)))-1.)
                        +8.1328E-03*(pow(10,-3.49149*((373.16/(T+273.16))-1))-1.)
                        +log10(1013.246));
        return(vp);
    }
    else
        return(-999.);
}

/*** vapour pressure over ice [mb] ******************************************/
double e_i_Calc(double T /* oC */)
{
   if (T > -273.16)
/* Goff-Gratch formula (List 1951) */
      return(pow(10,-9.09718*(273.16/(T+273.16)-1)
                    -3.56654*log10(273.16/(T+273.16))
                    +0.876793*(1-(T+273.16)/273.16)
                    +log10(6.1071)));
   else
      return(-999.);
}


/* The following equations are taken from Irabarne and Godson (1973)
   References are by chapter-equation                    */

/*** mixing ratio [] ********************************************************/
double r_Calc(double p    /* mb */,
              double T_d  /* oC */)
{
   if (T_d > -273.16 && p > 0)
      return(epsilon*e_w_Calc(T_d)/p);                /* IV - 77 */
   else
      return(-999.);
}

/*** dew point temperature [oC] *********************************************/
double T_d_Calc(double T  /* oC */, 
                double rh /* % */)
{
   if (rh > 0)
      return(1/(-4.25E-4*log10(rh/100.)+1/(T+273.16))-273.16);    /* VII - 6 */
   else
      return(-999.);
}

/*** wet-bulb temperature [oC] **********************************************/
double T_w_Calc(double p   /* mb */, 
                double T   /* oC */, 
                double T_d /* oC */)
{
double T_w, diff, delta;
int i;
   if (T > -273.16 && T_d > -273.16 && p > 0)
   {
      i=0;
      T_w = (T+T_d)/2;
      diff=T-T_w;
      do
      {
         i++;
         delta = l_v_Calc(T_w)*epsilon/c_p/p*(e_w_Calc(T_w) - e_w_Calc(T_d));
         diff = diff/2;
         if (T - delta > T_w)
            T_w = T_w + diff/2;
         else
            T_w = T_w - diff/2;
      }
      while (T_w+diff != T_w && i < 100);
      return(T_w);
   }
   else
      return(-999.);
}

/*** dry adiabat line *******************************************************/
double adiabat(double T,    /* oC */
               double theta /* K */)
{
   if (T > -273.16 && theta > 0)
      return(1000*pow((T+273.16)/(theta),1/chi_d));        /* VI - 1 */
   else
      return(-999.);
}

/*** Variation of T with height ********************************************/
double T_wrt_z(double T,  /* oC */
               double from_z, /* m */
               double to_z  /* m */)
{
   if (T > -273.16)
      return( T + 0.0065 * (from_z - to_z)); /* assume an SAO lapse rate */
   else
      return(-999.);
}

/*** Variation of rh with height ********************************************/
double rh_wrt_z(double T,  /* oC */
                double rh, /* % */
                double from_z, /* m */
                double to_z  /* m */)
{
double new_rh;

   if (T > -273.16 && rh > 0 && rh <= 100)
   {
      new_rh = 100 * e_w_Calc(T_d_Calc(T, rh))
             / e_w_Calc(T_wrt_z(T, from_z, to_z));
      if (new_rh > 100.) new_rh = 100.;   /* saturation reached */
      return(new_rh);
   }
   else
      return(-999.);
}

double rh_calc(double q, /* g/kg */
               double p, /* mb */
               double T  /* oC */)
{
double r, e, U;

    if (q > 0 && p > 0 && T > -273.16)
    {
        r = (q/1000.)/( 1.-(q/1000.) );    /* IV - 74 */
        e = p * r / (0.622 + r);        /* IV - 76 */

        if ( T > 0 )
            U = 100. * e/e_w_Calc(T);    /* IV - 83 */
        else
            U = 100. * e/e_i_Calc(T);    /* IV - 84 */

        return (U); 
    }
    else
      return(-999.);    
}
