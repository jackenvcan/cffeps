#ifndef DIURNAL_H
#define DIURNAL_H

typedef struct                        /* Climate data */
{
   double Tmax[12],
          Tmin[12],
          Tmean[12],
		  RH0600[12], 
		  W[12],
		  cured[12];

   int	  greenup[12], 
		  standing[12];

}  CLIM;

typedef struct                        /* Parameters from Beck & Trevitt 1989 */ 
{

   double alphaT,     /* hours from sunrise to T_min */
		  betaT,      /* hours from noon to T_max */
		  gammaT,     /* nocturnal decay factor */
		  alphaW,     /* hours from sunrise to W_min */
		  betaW,      /* hours from noon to W_max */
		  gammaW;     /* nocturnal decay factor */

}  BECK;


int ReadClimateData(CLIM * ptr);
int ReadBeckData(BECK * ptr);


		  
/************************ diurnal.h *****************************/

/*    These routines capture the various diurnal processes currently
      being used in fire modelling at the NoFC.                 */

/************************ hrffmc.c *****************************/

/* The official hourly FFMC routine based upon

   Van Wagner, C.E.  1977.  A method of computing fine fuel moisture 
      content thoughout the diurnal cycle.  Can. For. Serv.  Petawawa 
      For. Exp. Sta., Chalk R., Ont.  Inf. Rep. PS-X-38.

   Inputs include:
      temp      Temperature (oC)
      rh        Humidity (%)
      wind      Wind Speed (km/hr)
      rain      Last hour's rainfall (mm)
      oldffmc   Last hour's FFMC

   Return value is the new hour's FFMC.

*/
double hourly_ffmc(double temp, double rh, double wind, double rain, double oldffmc);

/************************ eqffmc.c *****************************/
/* An unofficial official hourly FFMC routine following the logic of

   Van Wagner, C.E.  1977.  A method of computing fine fuel moisture 
      content thoughout the diurnal cycle.  Can. For. Serv.  Petawawa 
      For. Exp. Sta., Chalk R., Ont.  Inf. Rep. PS-X-38.

   This routine calculates an FFMC in equilibrium with the environment.
   It is likely more closer to reality than CVW's original hourly calculations

   Inputs include:
      T    Temperature (oC)
      H    Humidity (%)
      W    Wind Speed (km/hr)
      rf   Last hour's rainfall (mm)
      F    Last hour's FFMC

   Return value is the new hour's FFMC.

*/
double eq_ffmc(double T, double H, double W, double rf, double F);

/************************ bdlffmc.c *****************************/

/* The function bdl_ffmc calculates the diurnal FFMC based upon

   Lawson, B.D., Armitage, O.B., Hoskins, W.D.  1996.  Dirunal Variation
      in the Fine Fuel Moisture Code: Tables and Computer Source Code.
      Canad-B.C. Partnership Agreement on For. Resour. Dev.: FRDA II, 
      Can. For. Serv., Victoria, British Columbia.  FRDA Rep. 245. 20 p.

   The function first checks that the program inputs are within the
   appropriate ranges. If any of the input values are out of
   range an error message is printed to the screen and control is
   returned to the operating system.

   NOTE: Appropriate ranges for inputs are defined as follows:
      hour    ( >= 1    and <= 2459  ) (integer)
      RH      ( >= 0    and <= 100   ) (integer)
      FF_FFMC ( >= 17.5 and <= 100.9 ) (one decimal place)
   Where: 
      FF_FFMC = 17.5  corresponds to the lower limit of Van Wagner's
	        original diurnal adjustment graphs;
      FF_FFMC = 100.9 corresponds to the theoretical upper limit of
                the current FF_FFMC scale.

    The function then calculates the moisture content at
    1600 (mc1600) based on the std. FF - FFMC value.
    This moisture content is then used in the appropriate equation
    to predict the moisture content at the desired time (adj_mc)
    between 1200 noon on the current day and 1159 the next day.
    Calculations for the various hours are based on the following:

       1200 to 2000 : Equations for every hour (interpolated for minute resolution).
       2001 to  559 : Interpolation between 2000 and 600 (using high RH equation)
       600  to 1100 : Equations for every hour (interpolated for minute resolution)
       1101 to 1159 : Extrapolation using 1100 and 1000 as end-points.

    The function then calculates the new FFMC from the new calculated
    moisture content and returns this time adjusted FFMC value.
*/
double bdl_ffmc(double ff_ffmc, int hour, int rh);

/************************ judi.c *****************************/

/*-------------------------------------------------------------------
   sunrise_calc was modified from a C program called sunrise.c 
   written by Bear Giles (bear@fsl.noaa.gov) with algorithms for solar
   declination, equation of time, and length of day provided by 
   Joseph Bartlo (jabartlo@delphi.com)

   NOTE: This routine is called by diurnal_calc and is included here 
   for completeness.  It does not need to be called any further to
   acheive J.A. Beck's diurnal curves.
   -------------------------------------------------------------------*/
int sunrise_calc (double lat,      /* Latitude (decimal degrees) */
                  double lng,      /* Longitude (decimal degrees) */
                  int    jday,     /* Julian Day */
                  double *sunrise, /* Local Standard Time (decimal hours) */
                  double *sunset   /* Local Standard Time (decimal hours) */);

/*-------------------------------------------------------------------
   This subroutine is based upon

      Beck, J.A.; Trevitt, C.F.  1989.  Forecasting diurnal variations in
         meteorological parameters for predicting fire behaviour.  
         Can. J. For. Res. 19: 791-797.

   equations numbers are included in remarks.

   Only deviation worth noting are:

1. Vapour pressure is passed as an input parameter.  It is assumed to be
   constant for the day.

2. Value and time of minimum temperature for the next day are assumed to 
   be same as today.

3. Equation 4 has been corrected (I think) to account for morning values.

   Inputs:
            lat     Latitude (decimal degrees)
            lng     Longitude (decimal degrees)
            jday    Julian Day
            t       Local Standard Time (decimal hours) 
            Tn      Minimum Temperature (Celsius)
            Tx      Maximum Temperature (Celsius)
            vp      Vapour Pressure (millibars)
            Wx      Maximum Wind Speed (kmh)
            Wn      Minimum Wind Speed (kmh)
            alphaT  Time lag between sunrise to time of min temp (hours)
            betaT   Time lag between solar noon to time of max temp (hours) 
            gammaT  A decay parameter from the paper (unitless)
            alphaW  Time lag between sunrise to time of min wind (hours)
            betaW   Time lag between solar noon to time of max wind (hours) 
            gammaW  A decay parameter from the paper (unitless)

   Outputs:
            T       Temperature at time t (Celsius)
            rh      Relative Humidity at time t (%)
            W       Wind Speed at time t (kmh)
  -------------------------------------------------------------------*/
void diurnal_calc(double lat, 
                  double lng, 
                  int jday, 
                  double t, 
                  double Tn, 
                  double Tx, 
                  double vp, 
                  double Wx, 
                  double Wn,
                  double alphaT, 
                  double betaT, 
                  double gammaT, 
                  double alphaW, 
                  double betaW, 
                  double gammaW, 
                  double *T, 
                  double *rh, 
                  double *W);

double EMCt(double t,
			double FFMCnoon, 
			double Tnoon, 
			double rhnoon, 
			double Wnoon, 
			double lat, 
			double lng, 
			int julian,
			double Tmin, 
			double Tmean, 
			double Tmax, 
			BECK *beck);

void Beck_calc(double lat, 
			   double lng, 
			   int jstart, 
			   int jend, 
			   double *Ts, 
			   double *Ws);

/************************ THE END *****************************/
#endif
