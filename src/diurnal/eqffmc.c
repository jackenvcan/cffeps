/*  EqFFMC.c

This routine calculates an FFMC in equilibrium with the environment.
It is likely more closer to reality than CVW's original hourly calculations

Equation numbers are based on:

Van Wagner, C.E., 1985: Equations and FORTRAN Program for the Canadian Forest Fire Weather Index System.
    Can. For. Serv., Ottawa, Ont. For. Tech. Rep. 33. 18 pp.

*/

/* modified Dec 23, 1996 for data type consistency */

#include <math.h>
#include <stdio.h>

double eq_ffmc(double T, double H, double W, double rf, double Fo)
{

double mo, m, mr, Ed, Ew, F;

    if (Fo <  0 || Fo > 101 ||
         T < -50.|| T > 50 ||
        H < 0   || H > 100 ||
        W < 0   || W > 200 ||
        rf < 0  || rf > 200)
        return (-99.);
    else
    {
        mo = 147.2*(101.-Fo)/(59.5 + Fo);                        /* 1 */

        Ed = 0.942 * pow(H,0.679) + 
             11. * exp((H-100.)/10.) + 
             0.18 * (21.1-T) * (1.-exp(-0.115*H));                /* 4 */

        Ew = 0.618 * pow(H,0.753) +
             10. * exp((H-100.)/10.) +
             0.18 * (21.1-T) * (1.-exp(-0.115*H));                /* 5 */

        m = mo;

/* Instantaneous recovery */
        if (mo > Ed)
            m = Ed;                                                /* 8 with the decay term removed */
        if (mo < Ew)
            m = Ew;                                                /* 9 with the decay term removed */


/* wetting phase is now conducted on equilibrium mc following FFMC calculations*/
        if (rf>0.)
        {
            mo = m;
            if (mo <= 150.) 
                mr = mo + 42.5*rf*exp(-100./(251.-mo))*(1.-exp(-6.93/rf));        /* 3a */
            else
                mr = mo + 42.5*rf*exp(-100./(251.-mo))*(1.-exp(-6.93/rf)) 
                        + pow(0.0015*(mo-150.), 2)*sqrt(rf);                    /* 3b */

            if (mr>250.) 
                mr=250.;

            m = mr;
        }

/* new FFMC */
        F = 59.5*(250.-m)/(147.2+m);                            /* 10 */

        if (F<0) printf("%lf %lf\n", F, m);
        return ( F );
    }
}
