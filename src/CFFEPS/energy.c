#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fbp_2009.h"
#include "feps_2011.h"

FILE *infile;

/* Energy Balance of a Fire required for Plume Rise

Anderson, K.R.; Pankratz, A; Mooney, C. 2011.  A thermodynamic approach to estimating smoke plume heights.
    In 9th Symp. on Fire and Forest Meteorology, Oct 18-20, 2011.  Palm Springs, CA.
        Am. Meteorol. Soc., Boston, MS.

Bulk densities for FBP fuel types are based on

Anderson, K.R.  2000.  Incorporating smoldering into fire growth modelling.
    Pages 31-36 in 3rd Symp. on Fire and Forest Meteorology, Jan. 9-14, 2000,
    Long Beach, CA.  Am. Meteorol. Soc., Boston, MS.

*/

double EnergyCalc(FEPS *feps, FBP *fbp)

/* inputs */
// int    header,            /* print header info 0=Yes, else No */
//          ddate,           /* detection date of fire YYYYMMDD */
//        dtime,           /* detection time of fire HHMM */
//          LDT,               /* local daylight time: 1=Yes, 0=No */
//          diurnal,           /* Diurnal adjustment: 0=Yes, else No */
//          print;           /* print: 1=Yes, 0=No */
//
// double obs,             /* observation time (time all values are assumed to be collected) */
//          temp,            /* surface temperature */
//          rh,              /* humidity */
//          DMC,
//          DC,
//          pressure,        /* MSL pressure */
//        growth,          /* change in fire size over last timestep */
//        area,            /* size of fire at dtime */
//        perimeter,       /* perimeter at dtime */
//          lapse,           /* lapse rate oC/km */
//          ts,              /* temperature at surface */
//          t850,            /* temperature at 850 mb */
//          t700,            /* temperature at 700 mb */
//          t500,            /* temperature at 500 mb */
//          zs,              /* height at surface */
//          z850,            /* height at 850 mb */
//          z700,            /* height at 700 mb */
//          z500,            /* height at 500 mb */
//          alpha,           /* entrainment half-angle */
//        Qo,              /* amount of energy previously injected into atmosphere */
//          timestep;           /* timestep for black boy radiation heat loss decimal hours */

{
int i, iFuel;
double     cpw = 4.1855 /* J g-1 K-1 */,
        cpwood = 1.7 /* J g-1 K-1 */,
        lv = 2.501E6 /* J kg^-1 */,
        H = 18000. /* heat of combustion [kJ kg^-1] */,
        pi = 3.1415926,
        sigma = 5.67E-8, /* Stefan Boltzmann Constant W m^-2 K^-4 8? */
        As, Aw, Ag, rs,  Mw, w, Qfire, Qw, Qf, Qs, Qr, Qinc, Fi, Qplume, L, mc, dT,
        Tfire, Ts, ps, A, h,
        rho2, rho4, rho6, rho8, DOB;

/* 1. Units are converted to SI units. */
        H = H*1000.; /* convert kJ kg^-1 to J kg^-1 */
        cpw = cpw*1000.; /* convert J g^-1 K^-1 to J kg^-1 K^-1*/
        cpwood = cpwood*1000.; /* convert J g^-1 K^-1 to J kg^-1 K^-1*/
        Tfire = 273.16 + + 800 + 400*fbp->CFB; /* [K] -- http://wildfiretoday.com/2011/02/26/at-what-temperature-does-a-forest-fire-burn/ */
        Ts = 273.16 + feps->temp; /* convert Ts to [K] */
        ps = feps->pressure * 100.;    /* convert ps to [pa] */
        As = A = 100. * 100. * feps->area; /* convert A to m^2 */
        Ag = 100. * 100. * feps->growth; /* convert area growth to m^2 */

/* 2. Determine fire area, growth, perimeter and radius. */

/* This assumes the entire fire perimeter is generating smoke */
    if (feps->perimeter > 0.)
        L = feps->perimeter;    /* assume a circular fire */
    else
        L = 2.*pi*sqrt(feps->area/pi);    /* assume a circular fire */

//    if (Ag <= 0.) Ag = A;   /* use the fire area in the case where growth is not used */
    if (Ag < 0.) Ag = A;   /* use the fire area in the case where growth is not used */

//  WRONG: does not apply here: if (feps->residencetime > 0. && feps->residencegrowth > 0. ) Ag = As = 100. * 100. * feps->residencegrowth; // new idea - this converts growth to and hourly growth and uses that to denote the fire size for plume volume calculation

    rs = sqrt(As/pi);            /* radius - eqn 18 (recalculated for case when Qplume was known */

/* 3. Assign default duff characteristics by fuel type based on the literature.  */

   for (i = 0; i < MAX_FUELS; i++){
      if (!strncmp(fbp->FuelType, Fuels[i], 2)) {
         iFuel= i;
         break;
      } else {
         iFuel = FUELTYPE_NF;
      }
    }
   switch (iFuel){
/* Use values from Anderson 2000 in g/cm^3 */
      case (FUELTYPE_C2):
         rho2 = 0.019; rho4 = 0.034; rho6 = 0.051; rho8 = 0.056;
         break;
      case (FUELTYPE_C3):
         rho2 = 0.015; rho4 = 0.020; rho6 = 0.032; rho8 = 0.066;
         break;
      case (FUELTYPE_C4):
         rho2 = 0.022; rho4 = 0.029; rho6 = 0.045; rho8 = 0.059;
         break;
      case (FUELTYPE_C6):
         rho2 = 0.030; rho4 = 0.050; rho6 = 0.050; rho8 = 0.050;
         break;
      case (FUELTYPE_C7):
         rho2 = 0.100; rho4 = 0.100; rho6 = 0.050; rho8 = 0.000;   /* values used to match SFC limits in FBP */
         break;
      case (FUELTYPE_M1):
         rho2 = (0.019/2); rho4 = (0.034+0.108)/2; rho6 = (0.051+0.108)/2; rho8 = (0.056+0.108)/2;  /* average of C2 and D1*/
         break;
      case (FUELTYPE_M2):
         rho2 = (0.019/2); rho4 = (0.034+0.108)/2; rho6 = (0.051+0.108)/2; rho8 = (0.056+0.108)/2;  /* average of C2 and D1*/
         break;
      case (FUELTYPE_M3):
         rho2 = 0.041; rho4 = 0.061; rho6 = 0.084; rho8 = 0.112;
         break;
      case (FUELTYPE_M4):
         rho2 = 0.041; rho4 = 0.061; rho6 = 0.084; rho8 = 0.112;
         break;
      case (FUELTYPE_S1):
         rho2 = 0.000; rho4 = 0.2; rho6 = 0.2; rho8 = 0.;  /* values used to match SFC limits in FBP */
         break;
      case (FUELTYPE_S2):
         rho2 = 0.000; rho4 = 0.5; rho6 = 0.3; rho8 = 0.;  /* values used to match SFC limits in FBP */
         break;
      case (FUELTYPE_S3):
         rho2 = 0.000; rho4 = 0.6; rho6 = 1.0; rho8 = 0.;  /* values used to match SFC limits in FBP */
         break;
      default:
       rho2 = rhoBs[iFuel];
       rho4 = rhoBs[iFuel];
       rho6 = rhoBs[iFuel];
       rho8 = rhoBs[iFuel];
   }

/* 4. Calculate depth of burn and the duff consumption from the surface fuel consumption (SFC) and densities */
   if (iFuel == FUELTYPE_01A || iFuel == FUELTYPE_01B) {
        mc = 147.2*(101.-fbp->FFMC)/(59.5+fbp->FFMC);
        Mw = mc/100. * fbp->GFL;
        DOB = 0.;
    }
    else
    {
        if (fbp->SFC < 2.*rho2*10.) /* rho*10 give mass per square metre per cm of depth; rho*20. gives mass per 2cm layer */
            DOB = fbp->SFC/rho2/10.;
        else if (fbp->SFC < rho2*20. + rho4*20.)
            DOB = 2 + (fbp->SFC - rho2*20.)/rho4/10.;
        else if (fbp->SFC < rho2*20. + rho4*20. + rho6*20.)
            DOB = 4 + (fbp->SFC - rho2*20-rho4*20)/rho6/10.;
        else
            DOB = 6 + (fbp->SFC - rho2*20-rho4*20-rho6*20.)/rho8/10.;


/* 5. Water mass is calculated from the moisture content as described by the FFMC, (top 2 cm), DMC (2-5 cm), DC (5+ cm) and the depth of burn.
      Add the water mass from the crown fuel consumption (TFC-SFC) and the foliar moisture content (FMC). */
        mc = 147.2*(101.-fbp->FFMC)/(59.5+fbp->FFMC);   // mc = moisture cotent

        if (DOB < 2.0)
            Mw = mc/100. * DOB/2.*(rho2*20.);           // Mw = mass of water
        else
        {
            Mw = mc/100. * rho2*20.;

            mc = exp((feps->DMC - 244.72)/(-43.43))+20.;

        if (DOB < 4.0)
                Mw = Mw + mc/100.*(DOB-2.)*(rho4*10.);
            else
            {
                Mw = Mw + mc/100.*(rho4*20.);
                if (DOB < 6.0)
                    Mw = Mw + mc/100.*(DOB-4.)*(rho6*10.);
                else
                {
                    Mw = Mw + mc/100.*(rho6*20.);
                    mc = 800./exp(feps->DC/400.);
                    Mw = Mw + mc/100.*(DOB-6.)*(rho8*10.);
                }
            }
        }
    }

    if (fbp->TFC > fbp->SFC)
        Mw = Mw + fbp->FMC/100. * (fbp->TFC - fbp->SFC);

    rs = sqrt(As/pi);            /* radius - eqn 18 (required for calculation of the fire wall) */
                                /* as written, this is independent of area growth */

    w = fbp->TFC;

/********************************************************************/

/* 6. Calculate Qfire, Qw, Qf, Qs, Qr, Qinc and from these Qplume  */

    Qfire = H * w * Ag;            /* eqn 2 - based on area growth */

    dT = 100. - feps->temp;
    Qw = Mw * (lv + cpw*dT) * Ag;    /* energy that goes into heating fuel moisture and evaporating it */

//    dT = (600 - 273.16) - feps->temp;    /* new value for heating temperature for fuel (~600 K) based on Albini et al 1995... but this is just the ignition temperature! */
    dT = 551.92 - feps->temp;    /* newer value for heating temperature for fuel (mid point between 600 K and 777oC) also based on Albini et al 1995 */

    Qf = w * cpwood * dT * Ag;    /* energy that goes into heating fuel to 500oC */

    Qs = 0.5 * H * fbp->SFC * Ag;  /* new approach assuming that 50% amount of heat equal to surface fuel consumption is injected into surface */

    Fi = fbp->CFB/2.;        /* assumption that the percent incomplete combustion is half of the Crown Fraction Burned */

//    feps->radiation = 2;

/* Radiation method 1 (default scheme) =- from Byram 1959, 1973: page 161  - assumption that 14% (17%?) lost due to radiation */
    Qr  = Qfire * 1200./8600.;

/* Radiation method 2 - Flame height to size of fire */
    if ((feps->radiation == 0) ||
        (feps->radiation == 2)) // Method 2
    {
        h = sqrt(fbp->HFI/300.);    /* Byram (1959, 1973 pg 175) */
        Aw = 2*pi*rs*h;             /* area of the wall surrounding the fire */

        if (Aw < As && As > 0.)
            Qr = Qfire * Aw/As/2;
    }

    if (feps->radiation == 3) // Method 3
//        Qr  = Qf + Qw; // new thinking... back radiation equal the forward radiation (assuming a steady state fire)
    {
        h = sqrt(fbp->HFI/300.);    /* Byram (1959, 1973 pg 175) */
        Tfire = 273.14 + 800 + 400*fbp->CFB; /* http://wildfiretoday.com/2011/02/26/at-what-temperature-does-a-forest-fire-burn/ */
        Qr = h * L * sigma * pow(Tfire, 4) * feps->timestep * 3600.;  /* assumption Qr is the backward radiation, which is based on the flame temperature */
    }

/* Radiation method 4 - use the value as a percentage of energy of the fire that enters the plume */
    if (feps->radiation > 3)
//        Qplume= feps->radiation /100. *Qfire;
        Qplume= feps->radiation /100. *Qfire;
    else
    {

        if (feps->print == 1)  printf(" %.2lf %.2e %.2e %.2e %.2e %.2e", feps->elapsed, Qfire, Qw, Qf, Qr, Qs);

/* Qinc is a weak link in the energy budget model */
//        Qinc = Fi * (Qf + Qw);    // new approach for Qinc - 10/08/2015 (Stull's eqn occasionally produced neg values)
        Qinc = Fi*Qf; // newer approach for Qinc - only energy of fuel burned is involved in incomplete combustion - 04-06-2018

        if (feps->print == 1)
            if (Qfire > 0.)
                printf(" %.2e ", Qinc );
//                printf(" %.2e ", Fi*(Qfire - Qw - Qf - Qr - Qs ));    /* Qinc */
            else
                printf(" %.2e ", Qfire);

/* relevant only for methods 1 and 2 */
        if (feps->sinks > 0)
//        Qplume = (1.-Fi)*(Qfire - Qw - Qf - Qr - Qs);   // old equation
//        Qplume = Qfire - (1 + Fi) * (Qw + Qf) - Qr - Qs;    // corrected equation based on discussions with Roland Stull and Rosie Howard - 05/28/2015
        Qplume = Qfire - Qw - Qf - Qr - Qs - Qinc;    // corrected equation based on discussions with Roland Stull and Rosie Howard - 05/28/2015
        else
            Qplume = Qfire;
    }
    if (Qplume<=0) Qplume=-999.;
    if (Qplume<=0) Qplume=0;

    feps->Qfire = Qfire;
    feps->Qplume = Qplume;

// if (Qplume*Qfire>0) printf ("... Qplume/Qfire = %lf\n",  Qplume/Qfire);

// printf("%.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e\n", Qfire, Qw, Qf, Qr, Qs, Qinc, Qplume);
//    printf("%.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, ", Qfire, Qw, Qf, Qr, Qs, Qinc, Qplume);

    return (Qplume);

}
