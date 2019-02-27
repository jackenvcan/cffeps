#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fbp_2009.h"
#include "feps_2011.h"

/* Emissions rates are based on CONSUME 3.0

Page numbers refer to:

    Consume 3.0 User's Guide
    Susan J. Prichard, Roger D. Ottmar, and Gary K. Anderson
    2004?
    
Bulk densities for FBP fuel types are based on
 
Anderson, K.R.  2000.  Incorporating smoldering into fire growth modelling.  
    Pages 31-36 in 3rd Symp. on Fire and Forest Meteorology, Jan. 9-14, 2000, 
    Long Beach, CA.  Am. Meteorol. Soc., Boston, MS.
    
*/

void EmissionsCalc(FEPS *feps, FBP *fbp)
{
double rho2, rho4, rho6, rho8, L, F, H, rhoL, rhoF, rhoH, LFC, FFC, HFC, DOB, DOB_F, DOB_S;
double Flaming, Smoldering, Residual;
int i, iFuel;

    Flaming = 0.;
    Smoldering = 0.;
    Residual = 0.;
    
    iFuel=-999;
   
/* Default duff characteristics by fuel type */
    for (i=0; i<17; i++)
        if (!strncasecmp(fbp->FuelType, Fuels[i], 2))  iFuel= i;

    if (iFuel>=0 && iFuel< 17 && feps->growth > 0. && fbp->TFC > 0.)
    {
        
/* 1. Apply the average bulk density per fuel type to bulk densities for 2, 4 6 and 8 cm depths */

/* these values are in fbp_2009.h */
        rho2 = rhoBs[iFuel];
        rho4 = rhoBs[iFuel];
        rho6 = rhoBs[iFuel];
        rho8 = rhoBs[iFuel];

/* 2. Overwrite these values with depth specific values for certain fuels. */

/* Use values from Anderson 2000 in g/cm^3 */
        if (!strncasecmp(fbp->FuelType, "C1", 2)) { rho2 = 0.045; rho4 = 0.045; rho6 = 0.045; rho8 = 0.045; }
        if (!strncasecmp(fbp->FuelType, "C2", 2)) { rho2 = 0.019; rho4 = 0.034; rho6 = 0.051; rho8 = 0.056; }
        if (!strncasecmp(fbp->FuelType, "C3", 2)) { rho2 = 0.015; rho4 = 0.020; rho6 = 0.032; rho8 = 0.066; }
        if (!strncasecmp(fbp->FuelType, "C4", 2)) { rho2 = 0.022; rho4 = 0.029; rho6 = 0.045; rho8 = 0.059; }
        if (!strncasecmp(fbp->FuelType, "C5", 2)) { rho2 = 0.093; rho4 = 0.093; rho6 = 0.093; rho8 = 0.093; }        
        if (!strncasecmp(fbp->FuelType, "C6", 2)) { rho2 = 0.030; rho4 = 0.050; rho6 = 0.050; rho8 = 0.050; }
        if (!strncasecmp(fbp->FuelType, "C7", 2)) { rho2 = 0.100; rho4 = 0.100; rho6 = 0.050; rho8 = 0.000; } /* values used to match SFC limits in FBP */
        
        if (!strncasecmp(fbp->FuelType, "D1", 2)) { rho2 = 0.061; rho4 = 0.061; rho6 = 0.061; rho8 = 0.061; } /* in the CFFEPS documentation */

/* not sure where these came from but they are in the CCFFEPS documentation */
        if (!strncasecmp(fbp->FuelType, "M1", 2)) { rho2 = (0.019+0.034)/2; rho4 = (0.034+0.108)/2; rho6 = (0.051+0.108)/2; rho8 = (0.056+0.108)/2; } /* average of C2 and D1*/
        if (!strncasecmp(fbp->FuelType, "M2", 2)) { rho2 = (0.019+0.034)/2; rho4 = (0.034+0.108)/2; rho6 = (0.051+0.108)/2; rho8 = (0.056+0.108)/2; } /* average of C2 and D1*/

        if (!strncasecmp(fbp->FuelType, "M3", 2)) { rho2 = 0.041; rho4 = 0.061; rho6 = 0.084; rho8 = 0.112; }
        if (!strncasecmp(fbp->FuelType, "M4", 2)) { rho2 = 0.041; rho4 = 0.061; rho6 = 0.084; rho8 = 0.112; }
        
        if (!strncasecmp(fbp->FuelType, "S1", 2)) { rho2 = 0.000; rho4 = 0.2; rho6 = 0.2; rho8 = 0.; } /* values used to match SFC limits in FBP */
        if (!strncasecmp(fbp->FuelType, "S2", 2)) { rho2 = 0.000; rho4 = 0.5; rho6 = 0.3; rho8 = 0.; } /* values used to match SFC limits in FBP */
        if (!strncasecmp(fbp->FuelType, "S3", 2)) { rho2 = 0.000; rho4 = 0.6; rho6 = 1.0; rho8 = 0.; } /* values used to match SFC limits in FBP */
   
/* old numbers */
        L = 2.;
        F = 5. - L;
        H = 25 - F - L;
        rhoL = rho2;
        rhoF = (2*rho4+rho6)/3.;
        rhoH = (rho6+2*rho8)/3.;

/* 3. Calculate the bulk densities of the L (0-1.2 cm), F (1.2-7 cm) and H (7-18 cm) layers from the 2 cm depth data. */
        
/* Nominal Fuel Depths [cm] from Van Wagner 1987, Table 1. */
        L = 1.2;
        F = 7. - L;
        H = 18. - F - L;
        rhoL = rho2;
        rhoF = ((2.-L)*rho2 + 2*rho4 + 2*rho6 + rho8)/F;
        rhoH = rho8;

/* 3.a. If grass fuel type (O1), burn off the entire grass fuel load. Use the grass reduction values for flaming, smoldering and residual combustion (see Table 2). */
        if (!strncasecmp(fbp->FuelType, "O1", 2)) /* These values are based on a personal communication with Bill de Groot */    
        {
            Flaming = 0.95*fbp->TFC;
            Smoldering = 0.05*fbp->TFC;
            Residual = 0.;

            DOB = DOB_F = DOB_S = 0.;    /* not calculated, not needed */
            LFC = FFC = HFC = 0.;    /* not calculated, not needed */
        }

        else
        {
            
/* 3.b. else, if slash, use the total fuel consumption (TFC). Use the slash reduction values for flaming, smoldering and residual combustion (see Table 2). */

/* Slash  - p 194*/
            if (!strncasecmp(fbp->FuelType, "S", 1))
            {
                Flaming = 0.70*fbp->TFC;
                Smoldering = 0.15*fbp->TFC;
                Residual = 0.15*fbp->TFC;

                LFC = FFC = HFC = 0.;    /* not calculated, not needed */
            }

            else
            {

/* 3.c.    else, burn off the LFH layers in sequence until the entire SFC is accounted for.  Then use the crown fuel consumption (CFC) for any canopy burned.  Use the corresponding reduction values for flaming, smoldering and residual combustion (see Table 2). */

                if (fbp->SFC < L*rhoL*10.) /* rhoL*10 give mass per square metre per cm of depth; rho*20. gives mass per 2cm layer */
                    LFC = fbp->SFC;
                else
                    LFC = L*rhoL*10.;

                if (fbp->SFC-LFC < F*rhoF*10.) /* rho*10 give mass per square metre per cm of depth; rho*20. gives mass per 2cm layer */
                    FFC = fbp->SFC-LFC;
                else
                    FFC = F*rhoF*10.;
            
                if (FFC < 0.) 
                    FFC = 0.;

                HFC = fbp->SFC-LFC-FFC;

                if (HFC < 0.) 
                    HFC = 0.;
            
/* Ground Fuels  - p. 166 */
                Flaming    = 0.90*LFC + 0.10*FFC + 0.00*HFC;
                Smoldering = 0.10*LFC + 0.70*FFC + 0.20*HFC;
                Residual   = 0.00*LFC + 0.20*FFC + 0.80*HFC;
        

/* Trees  - p. 169 */
                Flaming    = Flaming   + 0.80*(fbp->TFC-fbp->SFC) / 0.85;   /* using mid story as compromise between over/mid and understory */
                Smoldering = Smoldering + 0.05*(fbp->TFC-fbp->SFC) / 0.85;  /* " / 0.85 " used to equate Flaming + Smouldering to CFC */
            }
/* 4. Calculate depth of burn [cm] per stage */
            DOB = fbp->SFC/(rhoL*10.);
            if (DOB > L) DOB = L + (fbp->SFC - L*rhoL*10.)/(rhoF*10.);
            if (DOB > L+F)  DOB = L + F + (fbp->SFC - L*rhoL*10. - F*rhoF*10.)/(rhoH*10.);

/* not quite but close enough */
            DOB_F = Flaming/(rhoL*10);
            if (DOB_F > L) DOB_F = L + (Flaming - L*rhoL*10.)/(rhoF*10.);
            if (DOB_F > L+F)  DOB_F = L + F + (Flaming - L*rhoL*10. - F*rhoF*10.)/(rhoH*10.);

            DOB_S = DOB-DOB_F;

        }
    
// Flaming/Smoldering/Residual amounts in kg/m^2
        feps->Flaming = Flaming;
        feps->Smoldering = Smoldering;
        feps->Residual = Residual;

        feps->DOB = DOB;
        feps->DOB_F = DOB_F;
        feps->DOB_S = DOB_S;

/* 5. Calculate energy release per stage */        
        feps->Qflaming = feps->Qplume * Flaming / fbp->TFC;
        feps->Qsmoldering = feps->Qplume * Smoldering / fbp->TFC;
        feps->Qresidual = feps->Qplume * Residual / fbp->TFC;
        
    }
    else
    {
        rhoL = rhoF = rhoH = 0.;    /* not calculated, not needed */
        LFC = FFC = HFC = 0.;    /* not calculated, not needed */
        
        feps->Flaming = 0.;
        feps->Smoldering = 0.;
        feps->Residual = 0.;

        feps->DOB = 0.;
        feps->DOB_F = 0.;
        feps->DOB_S = 0.;
        
        feps->Qflaming = 0.;
        feps->Qsmoldering = 0.;
        feps->Qresidual = 0.;
    }
}

/* new code used in CFFEPS.c - incorporates the residence time per stage */

/*  In this approach, the residence times in each of the three are passed as fractions of an hour through the feps structure.
    It assumes that each phase is completed before the next begins
    It also assumes that the emissions rate in each phase is constant during the phase */

int EmissionsOverTime(double *Fs, double *Ss, double *Rs, int jModel, int iSpecies, EMISSIONS *emissions, FEPS *feps )
{
   double Fr, Sr, Rr, Ft, St, Rt, Fp, Sp, Rp, n, a, b, c, t1, t2;
   int i, imax;
   
   Fr = Sr = Rr = 0.;

/* 1. Calculate the total residence time in hours for flaming, smoldering and residual stages */
//     recalculated from hours to timesteps 2018-07-30
    imax = (int)ceil( (emissions->residence_flaming + emissions->residence_smoldering + emissions->residence_residual) / feps->timestep);
    if (imax > MAX_TIMESTEPS) imax = MAX_TIMESTEPS;
//    printf("imax = %d %lf %lf %lf %lf\n", imax, feps->timestep, emissions->residence_flaming, emissions->residence_smoldering, emissions->residence_residual);
    
    if (jModel >= 0)
    {
        a = emissions->a[jModel][iSpecies];
        b = emissions->b[jModel][iSpecies];
        c = emissions->c[jModel][iSpecies];
   
/* 2. Calculate the emissions rates [tonnes/hr] (or energy if jModel <0) for flaming, smoldering and residual stages */

//   tonnes/hr =             ha/timestep  / hours/timestep *      tonnes/ha
        if (emissions->residence_flaming > 0.)
            Fr = a / 1000. * feps->growth / feps->timestep * 10.*feps->Flaming;// / emissions->residence_flaming;   // was a/2000 for [lbs/ton] - KRA 2017-02-06 
        else
            Fr =0.;
        
        if (emissions->residence_smoldering > 0.)
            Sr = b / 1000. * feps->growth / feps->timestep  * 10.*feps->Smoldering;// / emissions->residence_smoldering; // was b/2000
        else
            Sr = 0.;

        if (emissions->residence_residual > 0.)
            Rr = c / 1000. * feps->growth / feps->timestep  * 10.*feps->Residual;// / emissions->residence_residual; // was c/2000
        else
            Rr = 0.;    

/* emissions factors a, b and c are now assumed to be in [g/kg] (not [lbs/ton] as described below) - KRA 2017-02-06 */

/*                 [lbs/ton]     [lbs/ton]     [ha]         [kg/m^2]         [hours]        
                           
   a, b and c are in lbs/ton    (changed to [g/kg] - 2017-02-06)
   divide by 2000 to convert lbs/ton to unitless (cancels units of a, b, c)
   feps->growth is in hectares
   feps->Flaming is a fraction of the TFC in kg/m^2 (coverted to tonnes per hectare by multiplying by 10)
   emissions->residence_flaming[iSpecies]]

   answer is in [metric] tonnes per hour
*/
    }
    else    /* if jModel <0, use this technique to come up with hourly energy per hour */
    {
//               J/timesteps    / hrs/timestep 
        if (emissions->residence_flaming > 0.)
            Fr = feps->Qflaming / feps->timestep;// / emissions->residence_flaming;
        else
            Fr =0.;

        if (emissions->residence_smoldering > 0.)
            Sr = feps->Qsmoldering / feps->timestep;// / emissions->residence_smoldering;
        else
            Sr = 0.;

        if (emissions->residence_residual > 0.)
            Rr = feps->Qresidual / feps->timestep;// / emissions->residence_residual;    
        else
            Rr = 0.;
    }

/* Resident times (hours) */
    Ft = emissions->residence_flaming;
    St = emissions->residence_smoldering;
    Rt = emissions->residence_residual;

/* Partial Resident times (hours) */   
    Fp = Sp = Rp = 0.;

/* 3. Stepping through each timestep */
    if (Fr+Sr+Rr > 0)
        for (i=0; i<imax;i++)
        {
            t1 = i * feps->timestep;
            t2 = t1 + feps->timestep;
/* 3.a. calculate partial residence times remaining in decimal hours */
            if (t2 <= Ft)
                Fp = feps->timestep;
            else 
                if (t1 < Ft)
                    Fp = modf(Ft-t1, &n);
                else
                    Fp = 0.;

            if (t2 <= Ft + St)
                Sp = feps->timestep - Fp;
            else 
                if (t1 < Ft + St)
                    Sp = modf(Ft+St-t1, &n) - Fp;
                else
                    Sp = 0.;

            if (t2 <= Ft + St + Rt)
                Rp = feps->timestep - Sp - Fp;
            else 
                if (t1 < Ft + St + Rt)
                    Rp = modf(Ft+St+Rt-t1, &n) - Sp - Fp;
                else
                    Rp = 0.;

/* 3.b. calculate the emissions for the partial hour */
//     recalculated from hours to timesteps 2018-07-30 

// tonnes/timestep = tonnes/hr * hrs * hrs/timestep / hrs
            if (emissions->residence_flaming > 0.)
                Fs[i] = Fr*Fp * feps->timestep / emissions->residence_flaming;
            else
                Fs[i] = 0;            
        
            if(emissions->residence_smoldering > 0.)
                Ss[i] = Sr*Sp * feps->timestep / emissions->residence_smoldering;
            else
                Ss[i] = 0.;

            if (emissions->residence_residual > 0.)
                Rs[i] = Rr*Rp * feps->timestep / emissions->residence_residual;
            else
                Rs[i] = 0.;
        }
    
    return(imax);
                
}

