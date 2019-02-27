#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "diurnal.h"
#include "thermo.h"


/*-------------------------------------------------------------------
   This subroutine was modified from a C program called sunrise.c 
   written by Bear Giles (bear@fsl.noaa.gov) with algorithms for solar
   declination, equation of time, and length of day provided by 
   Joseph Bartlo (jabartlo@delphi.com)
   -------------------------------------------------------------------*/

int sunrise_calc (double lat, double lng, int jday, 
                  double *sunrise, double *sunset)
{
double dtr, x, decl, eot, noon, t, daylength;

   if (lat<-80. || lat>80. || lng<-180. || lng>180.) exit (-1);

/* degrees to radians */
   dtr = 3.1415926/180.;

/* fraction of a year */
   x = dtr * 360./366. * (float) jday;

/* solar declination */
   decl = dtr * (0.33029
        - 22.9717 * cos (x) + 3.8346 * sin (x)
        - 0.3495 * cos (2. * x) + 0.0261 * sin (2. * x)
        - 0.1392 * cos (3. * x) + 0.0727 * sin (3. * x));

/* equation of time ? I don't recognize the significance of this */
   eot = -0.0001062 + 0.009554 * cos (x) - 0.122671 * sin (x)
         -0.053350 * cos (2. * x) - 0.156121 * sin (2. * x);

   noon = 12. - (lng / 15.) - eot;

/* ... so we'll just use Solar noon for now */
   noon = 12.;

/* fraction of a day in sunlight */
   t = -tan (dtr * (lat)) * tan (decl)
     + -sin (dtr * 0.8) / (cos (dtr * (lat)) * cos (decl));

/* 24 hrs night: sunrise & sunset = noon */
   if (t < -1. ) t = -1.;

/* 24 hrs daylight: sunrise & sunset = midnight */
   if (t > 1.)   t = 1.;

/* daylength in hours */
   daylength =  2. * (12. / 3.1415926) * acos (t);

   *sunrise = noon - daylength/2.;
   *sunset = noon + daylength/2.;

  return(0);
}
/***************************************************/
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

   -------------------------------------------------------------------*/

void diurnal_calc(double lat, double lng, int jday, double t, 
                  double Tn, double Tx, double vp, double Wx, double Wn,
                  double alphaT, double betaT, double gammaT, 
                  double alphaW, double betaW, double gammaW, 
                  double *T, double *rh, double *W)

/* Inputs:
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
*/
{
double tr, ts, tn, tx, 
       trnext, tsnext, tnnext,
       trprev, tsprev, tnprev, 
       ft, gt, Ts, Ws, svp;

int jprev, jnext;

/* Values for Denver */
   if (alphaT*betaT*gammaT == 0)
      alphaT = 0.88; betaT = 1.86; gammaT = -2.20;

   if (alphaW*betaW == 0)
      alphaW = 1.00; betaW = 1.24; gammaW = -3.59;

   jprev=jday-1; if (jprev<1)   jprev=365;
   jnext=jday+1; if (jnext>365) jnext=1;

    sunrise_calc(lat, lng, jprev, &trprev, &tsprev);   /* previous day */
    sunrise_calc(lat, lng, jday, &tr, &ts);
    sunrise_calc(lat, lng, jnext, &trnext, &tsnext);   /* next day */

/* Correct input time to between 0 and 24 */
    while (t > 24.) t = t-24.;
    while (t <  0.) t = t+24.;

/**** Temperature ****/
    if (Tn > -99. && Tx > -99.)
    {
        tnprev = trprev + alphaT;                          /* 5 */

        tn = tr + alphaT;                                  /* 5 */
        tx = 12. + betaT;                                  /* 6 */

        tnnext = trnext + alphaT;                          /* 5 */

/* Before sunrise */
        if (t <= tn)
        {
            ft = (tsprev - tnprev)/(tx - tnprev);           /* 3 */
            Ts = Tn + (Tx - Tn) * sin (3.1415926/2*ft);     /* 1 */

            gt = (24. + t - tsprev)/(24. - tsprev + tn);    /* 4a */
            *T = Tn + (Ts - Tn) * exp(gammaT * gt);         /* 2 */
        }

/* Between sunrise and sunset */
        if (tn < t && t < ts)
        {
            ft = (t - tn)/(tx - tn);                       /* 3 */
            *T = Tn + (Tx - Tn) * sin (3.1415926/2*ft);     /* 1 */
        }

/* After sunset */
        if (ts <= t)
        {
            ft = (ts - tn)/(tx - tn);                       /* 3 */
            Ts = Tn + (Tx - Tn) * sin (3.1415926/2*ft);     /* 1 */

            gt = (t - ts)/(24. - ts + tnnext);              /* 4 */
            *T = Tn + (Ts - Tn) * exp(gammaT * gt);         /* 2 */
        }

/**** Humidity ****/
        if (vp > 0.)
        {

        /*    Calculate saturation vapour pressure in mb */
/*            svp = 6.108*exp(17.27* (*T)/((*T) + 237.3)); */      /* 10 */

        /* a better routine */
            svp = e_w_Calc(*T);

            if (vp < svp)
                *rh = 100.* vp/svp;                             /* 9 */
            else
                *rh = 100.;
        }
    }

/**** Wind ****/
    if (Tn > -99. && Tx > -99.)
    {
        tnprev = trprev + alphaW;                          /* 5 */

        tn = tr + alphaW;                                  /* 5 */
        tx = 12. + betaW;                                  /* 6 */

        tnnext = trnext + alphaW;                          /* 5 */

        if (gammaW >= 0.)   /* sine - sine model */
        {

/* Before sunrise */
            if (t <= tn)
            {
                gt = (t - tx)/(tn + tx);                        /* 12 */
                *W = Wx - (Wx - Wn) * sin (3.1415926/2*gt);     /* 11 */
            }

/* Between sunrise and sunset */
            if (tn < t && t < ts)
            {
                ft = (t - tn)/(tx - tn);                       /* 3 */
                *W = Wn + (Wx - Wn) * sin (3.1415926/2*ft);     /* 11 */
            }

/* After sunset */
            if (ts <= t)
            {
                gt = (t - tx)/(tnnext + tx);                    /* 12 */
                *W = Wx - (Wx - Wn) * sin (3.1415926/2*gt);     /* 11 */
            }
        }
        else              /* sine - exp model */
        {
/* Before sunrise */
            if (t <= tn)
            {
                ft = (tsprev - tnprev)/(tx - tnprev);           /* 3 */
                Ws = Wn + (Wx - Wn) * sin (3.1415926/2*ft);     /* 1 */

                gt = (24. + t - tsprev)/(24. - tsprev + tn);    /* 4a */
                *W = Wn + (Ws - Wn) * exp(gammaW * gt);         /* 2 */
            }

/* Between sunrise and sunset */
            if (tn < t && t < ts)
            {
                ft = (t - tn)/(tx - tn);                       /* 3 */
                *W = Wn + (Wx - Wn) * sin (3.1415926/2*ft);     /* 1 */
            }

/* After sunset */
            if (ts <= t)
            {
                ft = (ts - tn)/(tx - tn);                       /* 3 */
                Ws = Wn + (Wx - Wn) * sin (3.1415926/2*ft);     /* 1 */

                gt = (t - ts)/(24. - ts + tnnext);              /* 4 */
                *W = Wn + (Ws - Wn) * exp(gammaW * gt);         /* 2 */
            }
        }
    }
}

/* This routine calculates the diurnal trend from noon weather and climate norms */

double diurnal_from_climate (double t, double Tnoon, double rhnoon, double Wnoon, 
                             double lat, double lng, int julian,
                             double Tmin, double Tmean, double Tmax, 
                             double Wmin, double Wmean, double Wmax, 
                             double alphaT, double betaT, double gammaT,
                             double alphaW, double betaW, double gammaW)
{
double Td, vp, noon, T, rh, W, Wx, Wn, Tx, Tn, EMCt, zero=0.;
int i;
   EMCt = 0.;
   if (Tnoon>0 && rhnoon>0)
   {
/* Calculate vapour pressure */
         Td = T_d_Calc(Tnoon, rhnoon);
         vp = e_w_Calc(Td);

/* Calculate the monthly average max and min winds (potentially a lot of error introduced by this) */
         Wx = 1.2 * Wnoon;
         Wn = .5 * Wnoon;

/* Estimate Tn, Tx starting from assumption that Tnoon = Tmean[imonth] */
         Tn = Tmin + (Tnoon - Tmean);
         Tx = Tmax + (Tnoon - Tmean);

/* Adjust Tn and Tx so Tnoon(obs) = Tnoon(pred) */
         noon = 12.;
         diurnal_calc(lat, lng, julian, noon, 
                      Tn, Tx, vp, Wx, Wn,
                      alphaT, betaT, gammaT, 
                      alphaW, betaW, gammaW, 
                      &T, &rh, &W);

         Tn = Tn + (Tnoon - T);   /* note that we are no longer interested in T */
         Tx = Tx + (Tnoon - T);

         if (Wnoon>0)
         {
            Wn = Wn * Wnoon/W;
            Wx = Wx * Wnoon/W;
         }
         else
         {
            Wn=0.; Wx=0.;
         }

/* If estimated min temp is less than observed dew point, 
      set min temp to observed dew point and 
      converge on max temp using min temp and noon temp */
         i = 0;
         if (Tn < Td)
         {
            Tn = Td;
            do
            {
               diurnal_calc(lat, lng, julian, noon, 
                            Tn, Tx, vp, Wx, Wn,
                            alphaT, betaT, gammaT, 
                            alphaW, betaW, gammaW, 
                            &T, &rh, &W);
               Tx = Tx + (Tnoon - T);
               i++;
            } while (fabs(Tnoon - T) > .1 && i < 100);
         }

         if (i == 100)
            printf ("Error converging on Tmax\n");
         else
         {
            diurnal_calc(lat, lng, julian, t, 
                         Tn, Tx, vp, Wx, Wn,
                         alphaT, betaT, gammaT, 
                         alphaW, betaW, gammaW, 
                         &T, &rh, &W);
         }
  }
   return(EMCt);
}

void Beck_calc (double lat, double lng, int jstart, int jend, double *Ts, double *Ws)
{
int i, j, tx, tn, tnn, wx, wn, wnn,
    ntn, ntx, ntg, nwn, nwx, nwg;
double Tn, Tx, Tnn, Tss, Wn, Wx, Wnn, Wss, tr, ts, trnext, tsnext,
       alphat, betat, gammat, alphaw, betaw, gammaw;

    ntn=0; ntx=0; ntg=0; nwn=0; nwx=0; nwg=0;
    alphat=0.; betat=0.; gammat=0.;
    alphaw=0.; betaw=0.; gammaw=0.;

    for (j=0; j<jend-jstart; j++)
    {
        sunrise_calc(lat, lng, j+jstart, &tr, &ts);
        sunrise_calc(lat, lng, j+jstart+1, &trnext, &tsnext);

        tx=0; tn=0; wx=0; wn=0;    /* max-min times */
        Tn=999; Tx=-999; Wn=999; Wx=-999; /* max-min values */

        Tnn=999; Wnn=999; tnn=0; wnn=0;     /* next day's min values */

/* Temp and wind at sunset */
        Tss = (*(Ts+24*j+(int)ts+1) - *(Ts+24*j+(int)ts)) * (ts-(int)ts) + *(Ts+24*j+(int)ts);
        Wss = (*(Ws+24*j+(int)ts+1) - *(Ws+24*j+(int)ts)) * (ts-(int)ts) + *(Ws+24*j+(int)ts);

//        printf("%lf %lf\n", ts, (ts-(int)ts));
//        printf("%lf %lf %lf\n", *(Ts+24*j+(int)ts) , *(Ts+24*j+(int)ts+1), Tss );
//        printf("%lf %lf %lf\n", *(Ws+24*j+(int)ts) , *(Ws+24*j+(int)ts+1), Wss );

        for (i=0;i<36;i++)
        {
            if ( *(Ts+24*j+i) < Tn && i < 12)
            {
                Tn = *(Ts+24*j+i);
                tn = i;
            }

            if ( *(Ts+24*j+i) < Tnn && i > ts)
            {
                Tnn = *(Ts+24*j+i);
                tnn = i;
            }

            if ( *(Ts+24*j+i) > Tx)
            {
                Tx = *(Ts+24*j+i);
                tx = i;
            }

            if ( *(Ws+24*j+i) < Wn && i < 12)
            {
                Wn = *(Ws+24*j+i);
                wn = i;
            }

            if ( *(Ws+24*j+i) < Wnn && i > ts)
            {
                Wnn = *(Ws+24*j+i);
                wnn = i;
            }

            if ( *(Ws+24*j+i) > Wx)
            {
                Wx = *(Ws+24*j+i);
                wx = i;
            }
        } /* for (i=0;i<36;i++) */

        if (tn > tr && tn < 12) 
        {
            alphat = alphat + tn - tr;
            ntn++;
        }

        if (tx > 12 && tx < ts)
        {
            betat = betat + tx - 12;
            ntx++;
        }

        if (wn > tr && wn < 12) 
        {
            alphaw = alphaw + wn - tr;
            nwn++;
        }

        if (wx > 12 && wx < ts)
        {
            betaw = betaw + wx - 12;
            nwx++;
        }

        for (i=(int)ts+1;i<tnn;i++)
        {
            if ((*(Ts+24*j+i) > Tnn)  && (Tss > Tnn) && (i > ts) && (24-ts-tnn) !=0)
            {
                gammat = gammat + log((*(Ts+24*j+i)-Tnn)/(Tss - Tnn))/((i - ts)/(tnn - ts));
                ntg++;
            }

            if ((*(Ws+24*j+i) > Wnn) && (Wss > Wnn) && (i > ts) && (24-ts-wnn) !=0)
            {
                gammaw = gammaw + log((*(Ws+24*j+i)-Wnn)/(Wss - Wnn))/((i - ts)/(wnn - ts));
                nwg++;
            }

        } /* for (i=(int)ts+1;i<36;i++) */

    } /* for (j=0; j<jend-jstart; j++) */

    alphat = alphat/ntn;
    betat = betat/ntx;
    gammat = gammat/ntg;

    alphaw = alphaw/nwn;
    betaw = betaw/nwx;
    gammaw = gammaw/nwg;

    printf("%lf %lf %lf %lf %lf %lf ", alphat, betat, gammat, alphaw, betaw, gammaw);
    printf("%d %d %d %d %d %d ", ntn, ntx, ntg, nwn, nwx, nwg);

}

/* A routine to compare results with those presented in Beck & Trevitt * /
void main()
{
double lat, lng, t, Tn, Tx, vp, Wx, Wn,
       alphaT, betaT, gammaT, 
       alphaW,  betaW,  gammaW, 
       T,  rh,  Wse, Wss;
int jday;

FILE *outfile;

/ * Gouldburn * /
   lat= -34.75;
   lng= 150.;
   alphaT = 0.50; betaT = 2.05; gammaT = -2.45;

/ * Feb 7 * /
   jday=38;

   Tn = 14.;
   Tx = 30.;
   vp = 10.;
   Wx = 25.;
   Wn = 5.;

   outfile = fopen("out.txt", "w");

   for (t==0.; t<24.; t++)
   {
          
      alphaW = 0.29; betaW = 2.91; gammaW = -2.89;
      diurnal_calc(lat, lng, jday, t,
                Tn,  Tx,  vp,  Wx,  Wn,
                alphaT,  betaT,  gammaT, 
                alphaW,  betaW,  gammaW,
                &T, &rh, &Wse);

      alphaW = 0.69; betaW = 1.57; gammaW = 0.;
      diurnal_calc(lat, lng, jday, t,
                Tn,  Tx,  vp,  Wx,  Wn,
                alphaT,  betaT,  gammaT, 
                alphaW,  betaW,  gammaW,
                &T, &rh, &Wss);

      
      fprintf(outfile, "%5.2lf %5.2lf %5.2lf %5.2lf %5.2lf\n", t, T, rh, Wse, Wss);
   }
   fclose(outfile);
} */
