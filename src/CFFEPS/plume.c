#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fbp_2009.h"
#include "feps_2011.h"

double r_Calc();
double l_v_Calc();

FILE *infile;

/* Plume model based on 

Anderson, K.R.; Pankratz, A; Mooney, C. 2011.  A thermodynamic approach to estimating smoke plume heights.  
    In 9th Symp. on Fire and Forest Meteorology, Oct 18-20, 2011.  Palm Springs, CA. 
    Am. Meteorol. Soc., Boston, MS.
    
*/

double PlumeCalc(FEPS *feps)

/* inputs - feps structure */
//   int  header,             /* print header info 1=Yes, 0=No */
//          ddate,            /* detection date of fire YYYYMMDD */
//        dtime,            /* detection time of fire HHMM */
//          LDT,                /* local daylight time: 1=Yes, 0=No */
//          diurnal,            /* Diurnal adjustment: 0=Yes, else No */
//          print;            /* print: 1=Yes, 0=No */
//
//   double obs,            /* observation time (time all values are assumed to be collected) */
//          temp,             /* surface temperature [oC] */
//          rh,               /* humidity [%] */
//          DMC,              /* Duff Moisture Code */
//          DC,               /* Drought Code */
//          pressure,         /* MSL pressure [mb] */
//        growth,           /* change in fire size over last timestep [ha] */
//        area,             /* size of fire at dtime [ha] */
//        perimeter,        /* perimeter at dtime [km] */ 
//          lapse,            /* lapse rate [oC/m] */
//          ts,               /* temperature at surface [oC] */
//          t850,             /* temperature at 850 mb [oC] */
//          t700,             /* temperature at 700 mb [oC] */
//          t500,             /* temperature at 500 mb [oC] */
//          zs,               /* height at surface [m] */
//          z850,             /* height at 850 mb [m] */
//          z700,             /* height at 700 mb [m] */
//          z500,             /* height at 500 mb [m] */
//          alpha,            /* entrainment half-angle [o] */
//        Qo,               /* amount of energy previously injected into atmosphere [j] */
//          timestep;            /* timestep for black boy radiation heat loss [decimal hours] */

{
int i, ii, imax;
double     g = 9.80025    /* gravitational acceleration [m sec^-2] */,
        cp = 1005. /* J kg^-1 K^-1 */,
        Rd = 287.05 /* J kg^-1 K^-1 */,
        Ld = -9.8 / 1000. /* dry adiabatic lapse rate [oC m^-1] */,
        pi = 3.1415926,
        sigma = 5.67E-8, /* Stefan Boltzmann Constant W m^-2 K^-4 8? */
        As, Aw, At, Ag, rs, rt, rho, Vi,  Mo, Mi, Mw, w, Ebb, Qbb, Qfire, Qw, Qf, Qs, Qr, Fi, Qt, dz, q, M, Qplume, Qatm, E, L, pt, mc, dT,
        Le, Tm, Te, Tt, Ts, ps, A, alpha, F, d, P, ddz, Mt, Ms, rhot, rhos;

/* 1. If the environmental lapse rate is unstable (<= -9.8 oC/km), or if the area or Qplume is less than zero, 
      the plume height is set to -9999 (missing) and the subroutine is exited */

    if (feps->lapse < -0.0098 || feps->area < 0. || feps->Qplume <= 0.) /* remove bad data */
    {
        dz = -9999.;
        M = -9999.;
// printf("Error: %lf %lf %lf\n", feps->lapse, feps->area, feps->Qplume);
    }
    else
    {
        
/* 2. Units are converted to SI units. */
        if (feps->lapse < -.0090) feps->lapse = -.009; /* lapse rate limited to -9.0 oC/km */
        P = feps->perimeter * 1000.;        /* convert perimeter to [m] */
        Le = feps->lapse;                   /* Environmental lapse rate [oC m^-1] */
        Ts = 273.16 + feps->ts;           /* convert Ts to [K] */
        ps = feps->pressure * 100.;            /* convert ps to [pa] */
        As = A = 100. * 100. * feps->area;  /* convert A to [m^2] */
        Ag = 100. * 100. * feps->growth;    /* convert area growth to [m^2] */
        alpha = feps->alpha * pi/180.;        /* convert to radians */

/* 3. Qplume is the energy injected into the plume as calculated by EnergyCalc */
        Qplume = feps->Qplume;

/* I don't use Area growth (per timestep) for anything in this routine */        

        if ( ( !strncasecmp(feps->shape, "line", 3) || !strncasecmp(feps->shape, "cresent", 3) )
            && feps->perimeter > 0.)
            d = As/P; /* forward spread distance [m] */
        else
            d = 0.;
            
/* This assumes the entire fire perimeter is generating smoke */
        if (feps->perimeter > 0.)
            L = P;    /* assume a circular fire [m] */
        else
            L = 2.*pi*sqrt(As/pi);    /* assume a circular fire [m] */
            
/*****  Testing this section *****/

/*  This approach turns out to be wrong
    as the plume height varies with the time step.  */
        
//        As = Ag;

/*  By keeping the plume volume based on the cumulative fire size, 
    the answer is stable regadless of time step.
    
    Note that energy is still calculated by area growth over time */

// new idea - this limited the fire size to growth over the residence time for plume volume calculation
        if (feps->residencetime > 0. && feps->residencegrowth > 0. ) As = 100. * 100. * feps->residencegrowth; 
    
        rs = sqrt(As/pi);            /* radius - eqn 18 (recalculated for case when Qplume was known */
        
/********************************************************************/        

/* 4. If past energy Qo present, add it to Qplume to get Qt */
        if (feps->Qo > 0.) 
            Qt = Qplume + feps->Qo;   /* previously injected heat of fire */
        else
            Qt = Qplume;

//        if (feps->print == 1) 
//            printf(" %.2e, %.2e\n", Qplume, Qt);

        dz = 1000.;     // first guest at height
        ddz = 1000.;    // step used in the iteractive process

        imax=10000;
        i=0;
     
        if (Qt <= 0.)
        {
           dz = 0.;
           M = 0.;
        }
        else
        {
            do 
            {
                rt = rs + dz*tan(alpha);                                        /* eqn 19 */
                At = pi*rt*rt;

/* 5.a. Calculate the energy per unit mass needed to heat plume to dry adiabat */
                q = -.5 * cp * Ld * dz * log(1. + dz * (Le - Ld) / Ts);            /* eqn 8 */
                
/* 5.b. Calculate the column mass of the plume */
                if (Le != 0.)
                {
                    pt = ps * (pow(1. + Le * dz / Ts, -g / Le / Rd) );              /* eqn 14 */ 
                    Mt = (ps - pt)/g * At;                                          /* eqn 15 */
                    Ms = (ps - pt)/g * As;                                          /* eqn 15 */
                    M = ps / g * As * (1. - pow(1. + Le * dz / Ts, -g / Le / Rd));    /* eqn 16 */
                }
                else    // in the isothermal case, take the average using Le+.5 and Le-.5
                {
                    pt = ps * ( (pow(1. + Le * dz / Ts, -g / Le+.5 / Rd))+(pow(1. + Le * dz / Ts, -g / Le-.5 / Rd) )/2);
                    Mt = (ps - pt)/g * At;                                          /* eqn 15 */
                    Ms = (ps - pt)/g * As;                                          /* eqn 15 */
                    M = ps / g * As * (1. - (pow(1. + Le * dz / Ts, -g / Le+.5 / Rd)+pow(1. + Le * dz / Ts, -g / Le-.5 / Rd))/2);                
                }
                Mo = M;  // Mo is mass without entrainment

                Tt = Ts + Le*dz; // plume top temperature
                rhot = pt / Tt / Rd;    // ideal gas law
                rhos = ps / Ts / Rd;    // ideal gas law

/* 5.c. Calculate air density in plume */
                if (Mt + Ms > 0)
                    rho = (Mt*rhot + Ms*rhos)/(Mt + Ms);    // new mass weighted scheme to determine rho
                else
                    rho = 0.;
                
                if (dz>0. && As>0.)
                    rho = Mo/As/dz;    // original scheme (need to use this to agree with paper exercise in FFM9 paper)
                else
                    rho = 0.;    // and later on M will e zero

/* 5.d. Adjust mass for entrainment */
                if (alpha > 0.)  /* include entrainment */
                {
                    M = 1./3.*pi*rho*dz*(rt*rt + rs*rt + rs*rs);                /* eqn 20 */                                            
                }

                if (!strncasecmp(feps->shape, "line", 3) )
                    M = rho * (d*dz*L + L*dz*dz*tan(alpha));        /* this calculation is for the volume of a wedge - no allowance for curvature */
                    
                if (!strncasecmp(feps->shape, "cresent", 3) )
                    M = rho * (d*dz*L + L*dz*dz*tan(alpha));        /* use a wedge following the burning perimeter (estimated in CFFEPS) */

/* 5.e. Calculate total energy needed to heat plume to dry adiabat */
                Qatm = q * M;                                                    /* eqn 17 */
                
/* 5.f. If energy required to heat the air mass (Qatm) is greater than plume energy (Qt), 
        halve the height-step size and reduce the height, otherwise increase the height */
                if (Qatm > Qt) 
                {
                    ddz = 0.5*ddz; // overshot the top, reduce the height-step size
                    dz = dz - ddz;
                }
                else
                    dz = dz+ddz;

                i++;
/* 5.g. Repeat calculations until the height-step size is less than a metre */
            } while ( fabs(ddz) > 1. && i < imax);  /* expect this to converge to within less than a metre */
        }

/* 6. From the previous plume energy Qo, calculate the black body radiation loss Qbb. */
        if (feps->sinks > 0)
            if (feps->timestep <= 0.)
                Qbb = 0.;
            else
            {
                Tt = Ts + Le*dz; // plume top temperature
                Tm = Ts +(Le - Ld)*dz; // modified surface temperature
                Qbb = feps->timestep * 3600. * dz * L * sigma * (pow( (Tm+Tt)/2, 4) - pow( (Ts+Tt)/2, 4));
            }
        else
            Qbb = 0.;
        feps->Qo = Qt - Qbb;

        if (feps->print == 1)  printf(" %.2e %.0lf\n", Qbb, dz);


        if (feps->Qo < 0.) feps->Qo = 0.;
   
    }

    if (dz>feps->zt && feps->zt>0.0001) dz = feps->zt;  // setting an upper limit equal to the tropopause (for whatever reason, feps->zt>0 fails whenn feps->zt=0.000

    feps->dz = dz;
    feps->M = M;

// printf ("PlumeCalc: %lf %lf %lf %lf %.2e %.2e %lf\n", q, M, pt/100., Tt-273.16, Qatm, Qt, dz);
    return (dz);
}


double PlumeCalcDry(FEPS *feps, UA *ua)

/*     This routine allows for a piecewise integration of an upper air profile to calculate plume rise
    It has been tested under the ICAO standard atmosphere and is within 1% of the PlumeCalc routine. 
    
    With that said, there appears to be an error in the calculations that I cannot trace -- apparently in energy calculations (q) and specifically in the temperature calculations (T1, T2, T3, T4) as absolute temperature ratio ln(theta1/theta2) is handled consistently.  Mass (M) is correct.  The calculation does not precisely measure the area in the trapezoid and, as a result, the addition of two adjacent trapazoids does not equal the area of the larger trapezoid.  The error tends to propagate through the calculations resulting in higher plumes with more piecewise calculations.  The error is small (tens of meters) but is annoying -- KRA 2018-05-16 */

/* inputs  - feps structure */
//   int  header,             /* print header info 1=Yes, 0=No */
//          ddate,            /* detection date of fire YYYYMMDD */
//        dtime,            /* detection time of fire HHMM */
//          LDT,                /* local daylight time: used as 1 = one hour offset */
//          diurnal,            /* Diurnal adjustment: 0=Yes, else No */
//          print;            /* print: 1=Yes, 0=No */
//
//   double obs,            /* observation time (time all values are assumed to be collected) */
//          temp,             /* surface temperature [oC] */
//          rh,               /* humidity [%] */
//          DMC,              /* Duff Moisture Code */
//          DC,               /* Drought Code */
//          pressure,         /* MSL pressure [mb] */
//        growth,           /* change in fire size over last timestep [ha] */
//        area,             /* size of fire at dtime [ha] */
//        perimeter,        /* perimeter at dtime [km] */ 
//          lapse,            /* lapse rate [oC/m] */
//          ts,               /* temperature at surface [oC] */
//          t850,             /* temperature at 850 mb [oC] */
//          t700,             /* temperature at 700 mb [oC] */
//          t500,             /* temperature at 500 mb [oC] */
//          zs,               /* height at surface [m] */
//          z850,             /* height at 850 mb [m] */
//          z700,             /* height at 700 mb [m] */
//          z500,             /* height at 500 mb [m] */
//          alpha,            /* entrainment half-angle [o] */
//        Qo,               /* amount of energy previously injected into atmosphere [j] */
//          timestep;            /* timestep for black boy radiation heat loss [decimal hours] */

/* inputs upper air sounding (passed as pointers to arrays of values) */
//  double *uT  temperature [oC]
//  double *uZ  height [m]
//  double *uP  pressure [mb]
//  int iLevels total number of levels in the array
{
int i, ii, iConverge;
double//     g = 9.80025,        /* gravitational acceleration [m sec^-2] */
        g = 9.80665,        /* Standard Gravity Wikipedia value */
        cp = 1005.,          /* J kg^-1 K^-1 */
        Rd = 287.05,        /* J kg^-1 K^-1 */
        Re = 6371000.,
        Ld = -9.8 / 1000.,    /* dry adiabatic lapse rate [oC m^-1] */
//        Ld = -g / cp,        /* dry adiabatic lapse rate [oC m^-1] */
        pi = 3.1415926,
        sigma = 5.67E-8,     /* Stefan Boltzmann Constant W m^-2 K^-4 8? */
        As, Aw, At, Ag, rs, rt, rho, Vi, Mo, Mi, Mw, w, Ebb, Qbb, Qfire, Qw, Qf, Qs, Qr, Fi, Qt, dz, q, M, m, Qplume, Qatm, E, L, pt, mc, dT,
        Le, Tm, Te, Tt, Ts, ps, A, alpha, F, d, P, dp, qc, Mc,
        T1, T2, T3, T4, Z1, Z2, Z3, Z4, P1, P2, P3, P4, Mt, Ms, rhot, rhos, ddz,
        theta1, theta2, theta3, theta4,    thetas, 
        T5, theta5, P5, Z5,    qa, qb, qT, ZZ, // new strategy where qT = qa+qb
        r1, r2, M1, M2, rho1, rho2, A1, A2;
        
double *T0, *Z0, *P0;

/* 1. If the area or Qplume is less than zero, the plume height is set to -9999 (missing) and the subroutine is exited  */
    if (feps->area < 0. || feps->Qplume <= 0.) /* remove bad data */
    {
        dz = -9999.;
        M = -9999.;
    }
    else
    {
/* 2. Units are converted to SI units. */
        if (feps->lapse < -.0090) feps->lapse = -.009; /* lapse rate limited to -9.0 oC/km */
        P = feps->perimeter * 1000.;        /* convert perimeter to [m] */
        Le = feps->lapse;                   /* Environmental lapse rate [oC m^-1] */
        As = A = 100. * 100. * feps->area;  /* convert A to [m^2] */
        Ag = 100. * 100. * feps->growth;    /* convert area growth to [m^2] */
        alpha = feps->alpha * pi/180.;        /* convert to radians */
        ZZ = 1.0;                            /* the vertical step used in converging [m] */

/* 3. Qplume is the energy injected into the plume as calculated by EnergyCalc */
        Qplume = feps->Qplume;
        
/* I don't use Area growth (per timestep) for anything in this routine */     
                               
        if ( ( !strncasecmp(feps->shape, "line", 3) || !strncasecmp(feps->shape, "cresent", 3) )
            && feps->perimeter > 0.)
            d = As/P; /* forward spread distance [m] */
        else
            d = 0.;
            
        /* This assumes the entire fire perimeter is generating smoke */
        if (feps->perimeter > 0.)
            L = P;    /* assume a circular fire [m] */
        else
            L = 2.*pi*sqrt(As/pi);    /* assume a circular fire [m] */
            
/*****  Testing this section *****/

/*  This approach turns out to be wrong
    as the plume height varies with the time step.  */
        
//        As = Ag;

/*  By keeping the plume volume based on the cumulative fire size, 
    the answer is stable regardless of time step.
    
    Note that energy is still calculated by area growth over time */
    
// new idea - this limited the fire size to growth over the residence time for plume volume calculation
        if (feps->residencetime > 0. && feps->residencegrowth > 0. ) As = 100. * 100. * feps->residencegrowth; 
 
        rs = sqrt(As/pi);            /* radius - eqn 18 (recalculated for case when Qplume was known */

        /********************************************************************/        

/* 4. If past energy Qo present, add it to Qplume to get Qt */
        if (feps->Qo > 0.) 
            Qt = Qplume + feps->Qo;   /* previously injected heat of fire */
        else
            Qt = Qplume;

//        if (feps->print == 1) 
//            printf(" %.2e, %lf, %lf %lf\n\n", Qplume, feps->M, feps->area, Qplume/feps->M);

/* 5. Set array pointers to first upper air data point. */
        T0 = ua->T;        // pointer to first record
        Z0 = ua->Z;        // pointer to first record
        P0 = ua->P;        // pointer to first record

/* Point 3 is the surface value (fixed) */
        T3 = *T0;       // temperature at/near surface [oC]
        Z4 = Z3 = *Z0;  // height at/near surface [m] -- constant
        P4 = P3 = *P0;  // pressure at/near surface [mb] -- constant

        Ts = T3+273.16;     // convert to [K]
        ps = P3*100.;       // convert to [pa]
        thetas = Ts;

/* 6. Set cumulative energy per unit mass, cumulative mass and plume height to zero. */
        qc = q = 0.;         // cumulative energy
        Mc = M = 0.;            // cumulative mass
        dz = 0.;
        iConverge = 0;  // track whether top has been reached (Qstm>Qplume) and convergence is required
        dp = 0.; dT = 0.; dz = 0.;
        
        if (Qt <= 0.)
        {
           dz = 0.;
           M = 0.;
        }
        else

/* 7. Begin converging on plume height. */
        {
            i=0;
            do
            {

/* these calculations are done in [oC] and [mb] */
/* 7.a. If stepping upwards (Qatm < Qt), read in next level of data */ 
                if (iConverge == 0)
                {
/* 7.a.i. Level 2 is the upper air data for the current level */
                    T2 = ua->T[i]; 
                    Z2 = ua->Z[i];
                    P2 = ua->P[i];

/* 7.a.ii Level 1 is the upper air data for the next level */
                    T1 = ua->T[i+1];
                    Z1 = ua->Z[i+1];
                    P1 = ua->P[i+1];

/* correct for geopotential height (plume rise is calculated in gpm) */
//                    Z1 = Z1 * Re*Re/(Z1 + Re)/(Z1+ Re);
//                    Z2 = Z2 * Re*Re/(Z2 + Re)/(Z2+ Re);

                    if (Z1!=Z2)
                        Le = (T1-T2)/(Z1-Z2);   // environmental lapse rate in upper layer
                    else
                        Le = 0.;
                    
                }
                else 

/* 7.b. Else level 2 remains the same; reduce Point 1 by 1 m, adjust pressure and temperature accordingly */
                {
                    if (Z1-ZZ <= Z2)  // not sure why this would happen but it seems to
                        iConverge = -2;
                    else
                    {
                        P1 = P2 + (P1-P2) * (Z1-Z2-ZZ)/(Z1-Z2);
                        T1 = T2 + (T1-T2) * (Z1-Z2-ZZ)/(Z1-Z2);
                        Z1 = Z1 - ZZ;
                    }
                }

                dz = Z1-Z4;

                rt = rs+dz*tan(alpha);  /* eqn 19 */
                At = pi*rt*rt;
                
                r2 = r1;                /* now calculate by layer */
                A2 = A1;
                r1 = rs+dz*tan(alpha);
                A1 = pi*r1*r1;

/* 7.c.    Calculate energy required to heat atmosphere to dry adiabat.

   Solution is found by calculating energy contained in trapezoid defined by (T1, theta1), (T2, theta2), (T3, theta3), (T4, theta4)
   where the temperatures are at the next level (T1), the current level (T2), 
   while T3 and T4 are those when lowered adiabatically to the surface (T1 to T4, T2 to T3) */

/* old approach - Irabarne and Godson Chapter VII - eqn 108 * /
                theta4 = theta1 = (T1+273.16)*pow((1000./P1), 0.2854);
                theta3 = theta2 = (T2+273.16)*pow((1000./P2), 0.2854);
                T3 = theta3*pow((P3/1000.), 0.2854) - 273.16;
                T4 = theta4*pow((P4/1000.), 0.2854) - 273.16; */
                
/* alternate approach to be consistent with PlumeCalc (note that theta is wrt the surface layer, not 1000 mb) * /
                theta4 = theta1 = (T1+273.16) - Ld*(Z1-Z4);
                theta3 = theta2 = (T2+273.16) - Ld*(Z2-Z3); */

/* 7.c.i. Calculate T3 and T4 using dry adiabat */
                T3 = T2 - Ld*(Z2-Z3);
                T4 = T1 - Ld*(Z1-Z4);
                T5 = T2 - Le*(Z2-Z3);

/* 7.c.ii. Calculate potential temperatures */
                theta1 = theta4 = (T4+273.16);
                theta2 = theta3 = (T3+273.16);
                theta5 = (T5+273.16);

/* 7.c.iii. If stepping upward, qc and Mc are updated (when stepping down, qc and MC remain constant) */
                if (iConverge == 0) 
                {    
                    qc = q;    // old energy
                    Mc = Mo;  // Mo is mass without entrainment
                }

// old strategy, calculate area in trapezoid -- don't know why this is wrong
//              q = qc - .5 * cp * log((theta1)/(theta2)) * (T1+T2-T3-T4);  // use the baseline qc as T1, P1, Z1 is converged on 0
    
/* 7.c.iv. Calculate area in trapezoid
          (new strategy, subtracting triangles to calculate area in trapezoid) */
                qT = -.5 * cp * log((theta1)/(theta5)) * (T1-T4);
                qa = -.5 * cp * log((theta2)/(theta5)) * (T2-T3); 
                qb = qT - qa;    // resulting area of trapezoid

/* 7.c.v. Add new energy and mass to previous values */
                q = qc + qb;
                Mo = Mc  + 100.*(P2 - P1) / g * As;                            /* eqn 15 - Mo is column mass*/

/* 7.d. Calculate air density in plume. */
                if (dz>0. && As>0.)
                    rho = Mo/As/dz;    // original scheme (need to use this to agree with paper exercise in FFM9 paper)
                else
                    rho = 0.;    // and later on M will e zero                if 
                rho = Mo/As/dz;

/* 7.e. Adjust mass for entrainment. */
                if (dz == 0) 
                    M=0.;
                else
                {        
                    if (alpha > 0.)  /* include entrainment */
                        M = 1./3.*pi*dz*(rt*rt + rs*rt + rs*rs) * rho;                /* 20 */

                    if (!strncasecmp(feps->shape, "line", 3) )
                        M = rho * (d*dz*L + L*dz*dz*tan(alpha));        /* this calculation is for the volume of a wedge - no allowance for curvature */
                    
                    if (!strncasecmp(feps->shape, "crescent", 3) )
                        M = rho * (d*dz*L + L*dz*dz*tan(alpha));        /* use a wedge following the burning perimeter (estimated in CFFEPS) */
                }
                
                Qatm = q * M;                                                    /* eqn 17 */

/* 7.f. If energy required to heat the air mass (Qatm) is greater than plume energy (Qt), 
        begin calculations downward, reducing the top level by 1 metre
        until Qatm < Qt again, in which case the soluiton is reached. */
                if (iConverge == 0 && Qatm < Qt && i < ua->iLevels && i < MAX_LEVELS)
                {
//                    uT++; uZ++; uP++;
                    i++;
                }
                else // the first time it fails will be when Qatm > Qplume (because iConverge will equal 0)
//                    if (iConverge == 1 && Qatm < Qplume)
                    if (iConverge == 1 && Qatm < Qt)
                        iConverge = -1; // done
                    else
                        if (iConverge == 0)
                            iConverge = 1;

            } while (iConverge >= 0);
        }

// printf ("PlumeCalcDry: %lf %lf %lf %lf %.2e %.2e %lf\n", q, M, P1, T1, Qatm, Qt, dz);

/* 6. From the previous plume energy Qo, calculate the black body radiation loss Qbb. */
        if (feps->sinks > 0)
            if (feps->timestep <= 0.)
                Qbb = 0.;
            else
            {
                Tt = Ts + Le*dz; // plume top temperature
                Tm = Ts +(Le - Ld)*dz; // modified surface temperature
                Qbb = feps->timestep * 3600. * dz * L * sigma * (pow( (Tm+Tt)/2, 4) - pow( (Ts+Tt)/2, 4));
            }
        else
            Qbb = 0.;
        feps->Qo = Qt - Qbb;

//        if (feps->print > 0)  printf(" %.2e %.0lf\n", Qbb, dz);

        if (feps->Qo < 0.) feps->Qo = 0.;
   
    }
    feps->dz = dz;
    feps->M = M;

    return (dz);
}
