#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include "fbp_2009.h"
#include "fwi84.h"
#include "feps_2011.h"

FILE *infile, *outfile, *csvfile, *profile;

double PlumeCalc();
double PlumeCalcDry();
double EnergyCalc();
void EmissionsCalc();
int EmissionsOverTime();
double bdl_ffmc();
double BUIcalc();
void FBPtest();
int mtime();
double Readprofile();
int Readupperair();
int Reademissions();
int ReadEmissionsData();
int ReadEmissionsModels();
int ReadFEPSData();
int WriteFEPSData();

int main(int argc, char *argv[])
{
    char infilename[NB_CHAR_PATH], outfilename[NB_CHAR_PATH], profilename[NB_CHAR_PATH], profilepath[NB_CHAR_PATH], buffer[480], errortext[100];
    double Le, Ld, Ts, Tm, Tt, Ths, Thm, q, m, surface, Ediff, Tplume, Tenv, lapse, dz, dH;
    int i, j, nshort, nlong, status, firsttime, iyear, imonth, iday, hour, batchjob, iLevels,
        pyear, pmonth, pday, phour, pmin, pDj;
    double x, y, FFMCy, FFMCnoon, profile_date, timezone, detectionUTC, detectionDj, sum,
       elapsed, detectionLST, discovery, ignitionLST, discoverysize, discoverydist, timestep, current, observationLST,
       t1, t2, dt, A, AA, P, D, dist, area, perim, nextdist, nextfdist, nextbdist, growth, a,b, diameter, hourlygrowth, dailygrowth,
       Alag, Plag, Dlag, tlag, precip,
       DC, DMC, FFMC, FFMCt, hum, mc, Q, Qo, Qfire;
    double Qs8, Qs7, Qs5, Q87, Q85, Les8, Les7, Les5, Le87, Le85, dzs8, dzs7, dzs5, dz87, dz85, dzb8, dzbs;
    double dz8, dz7, dz5, l8, l7, l5, r_smoke, estarea, estarea2, dznew, Mdry, Mplume;    
    double g = 9.80025, cp = 1005, sbc = 5.67e-8, H = 18000000.;
    double pi = 3.1415926;

    char rep_date[24], rep_time[20], cmc_date[20], source[6], sensor[6], fuel[12];
    char *token;
    double lat, lon, fwi, ros, sfc, tfc, bfc, hfi;
    double plat, plon, parea, ptime, totalemissions, briggs;    // used with the CMC data
    double Qs[MAX_HOURS+1]; // we need to consider both end points UTC(t=0) and UTC(t=MAX_HOURS)
    double *Fs, *Ss, *Rs;
    int iii, iUTC, iMax, iFuel, nFuels, iModel, jModel, iSpecies, nSpecies, cmcdate, cmcDj, cmcDjo, cmcUTC, cmcUTCo, cmcLST, newfire;
    
    double *areas, *areas0, *ts, *ts0;
    
    FBP fbp, pfbp;
    FEPS feps, pfeps;
    EMISSIONS emissions;
    EMODELS emodels;
    UA ua;

static int jdays[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};   

/* this is an artifical diurnal trend for fire growth */
static double weight[] = 
                  { 0.015648395,
                    0.012470576,
                    0.010137905,
                    0.008443153,
                    0.004159725,
                    0.003641906,
                    0.004455177,
                    0.005716244,
                    0.017972534,
                    0.032938579,
                    0.043543306,
                    0.055667556,
                    0.065778522,
                    0.078207403,
                    0.090168462,
                    0.104189127,
                    0.101294339, // adjusted to make the total = 1
                    0.091901405,
                    0.072798822,
                    0.057406612,
                    0.044004006,
                    0.033647664,
                    0.025825035,
                    0.019983547 };

/* 1. Set up the run */
// if (feps->print==-2)printf("1.\n");

/* 1.a. Read the argument line for the name of the input file.  
    If this is missing, use the default infile.dat. */                    
// if (feps->print==-2)printf("1.a.\n");
    strcpy(infilename, "infile.dat\0");

    if (argc > 1)
    {
        sprintf(infilename, "%s\0", argv[1]);
        strcpy(infilename, argv[1]);
    }

    
/* 1.b. setup FBP default values */
// if (feps->print==-2)printf("1.b.\n");
    FBPdefaults(&fbp); 
    
/* 1.c. read in FEPS input data */
// if (feps->print==-2)printf("1.c.\n");
    if (ReadFEPSData(infilename, &feps, &fbp)) printf("Errror reading %s", infilename);
    feps.perimeter = diameter = 2*sqrt(10000.*feps.area/2/pi)/1000.;    // diameter of a circle [km] of estarea [ha]
    estarea = feps.area;
    
    feps.radiation = 0;

/* 1.d. Define time series arrays for emissions of each pollutant */
// if (feps->print==-2)printf("1.d.\n");
    nSpecies = ReadEmissionsData(feps.emissionsfilename, &emissions);
    nFuels = ReadEmissionsModels(feps.emissionsmodelsfilename, &emodels);

    if (!strcasecmp(feps.method, "cmc") ) feps.timestep = 1.;
    
    lapse = feps.lapse; /* used to refresh the original lapse rate per time step */

    areas = (double *)calloc((int)(24./feps.timestep + 1), sizeof(double));  // used to track the areas 
    areas0 = areas;

    ts = (double *)calloc((int)(24./feps.timestep + 1), sizeof(double));  // used to track the areas 
    ts0 = ts;
    
    for (i=0; i<(int)(24./feps.timestep + 1); i++);
    {
        *areas = -999; areas++;
        *ts = -999; ts++;
    }
    areas = areas0; ts = ts0;
    

/* number of short/long-term smoldering time steps */
    nshort = (int)(2./feps.timestep);
    nlong = (int)(MAX_HOURS/feps.timestep);
      

    Fs = (double *)calloc((int)(nlong+1), sizeof(double));
    Ss = (double *)calloc((int)(nlong+1), sizeof(double));
    Rs = (double *)calloc((int)(nlong+1), sizeof(double));

    newfire = 0; // false
    profile_date = -999.;    
    totalemissions = 0.;
    dznew = 0.; // seems this is turned off for now

/*******************************************************************************************/
    if (!strcasecmp(feps.method, "cmc") )  /* CMC Firework method */
    {

/* 1.e. Open CMC GEM forecast CSV file as input.  Read in header line */
// if (feps->print==-2)printf("1.e.\n");
        if ( (csvfile = fopen(feps.csvfilename, "r")) == NULL) {
            printf("Error: could not find %s\n\n", feps.csvfilename);
        
        } else if ( (outfile = fopen(feps.outfilename, "w")) == NULL) {
            printf("Error: could not write to %s\n\n", feps.outfilename);        
        } 
        else 
        {
// p is shorthand for previous record 
            plat=0.;
            plon=0.;
            parea=0.;
            Qo=0.;
            
// read (and ignore) header info
            fgets(buffer, 480, csvfile);
            
// print header info
            if (!strcasecmp(feps.method, "cmc") ) // seems we're already in (!strcasecmp(feps.method, "cmc")
            {
                feps.timestep=1.;
                
/* 1.f. Print (to the screen), new header line for the CFFEPS predictions CSV file. */
// if (feps->print==-2)printf("1.f.\n");
                fprintf(outfile, "lat, lon, rep_date, source, sensor, ffmc, dmc, dc, ws, fwi, fuel, ros, sfc, tfc, bfc, hfi, estarea, ");
                fprintf(outfile, "UTC, temp, rh, ws, precip, TS, T850, T700, T500, T250, ZS, Z850, Z700, Z500, Z250, Area(t), Growth (t), TotalEmissions, r_smoke, Zplume, Mplume, Qo, Qplume");
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
                {
                    fprintf(outfile, ", %s F", emissions.species[iSpecies]);
                    fprintf(outfile, ", %s S", emissions.species[iSpecies]);
                    fprintf(outfile, ", %s R", emissions.species[iSpecies]);
                }
                fprintf(outfile, "\n");
            }

// read in first line of data            
            fgets(buffer, 480, csvfile);
        }  // if ( (csvfile = fopen(feps.csvfilename, "r")) == NULL)... else
        
// for "top-hat", reset weight[] values to top-hat approach
        if (!strncmp(feps.shape, "top-hat", 3))        
            for (i=0;i<24;i++)
               if (i>9 && i <=21)
                   weight[i] = 1./12.;
               else
                   weight[i] = 0.;

/* 2. Begin reading CFFEPS predictions CSV file */               
// if (feps->print==-2)printf("2.\n");
        do  /* while !feof(csvfile); */
        {
            
/* 2.a. Reset arrays to zero. */
// if (feps->print==-2)printf("2.a.\n");
        
// ... but I do this again in 2.d. ?
            for (i=0;i<nlong+1;i++)
            {
                Fs[i]=0.;
                Ss[i]=0.;
                Rs[i]=0.;
            }

            feps.header=-1;
            iii++;
            
/* 2.b. Read line of data into buffer.  Parse comma separated variables in buffer.  
        Make adjustments for possible variation in units (e.g Kelvin instead of Celsius, etc.) */
// if (feps->print==-2)printf("2.b.\n");

// lat, lon, rep_date, source, sensor, ffmc, dmc, dc, ws, fwi, fuel, ros, sfc, tfc, bfc, hfi, ...
            token = strtok(buffer, ",");
            fbp.LAT= atof(token);  token=strtok(NULL, ",");
            fbp.LON = atof(token);  token=strtok(NULL, ",");
            strcpy(rep_date, token);  token=strtok(NULL, ",");

if (feps.print<0 && strlen(rep_date)<16) printf("2.b. Warning: incomplete rep_date (%s),\n", rep_date);

            strcpy(source, token);  token=strtok(NULL, ",");
            strcpy(sensor, token);  token=strtok(NULL, ",");
            fbp.FFMC = atof(token);  token=strtok(NULL, ",");
            feps.DMC = atof(token);  token=strtok(NULL, ",");
            feps.DC = atof(token);  token=strtok(NULL, ",");
            fbp.WS = atof(token);  token=strtok(NULL, ","); /* overwritten later with CMC data */
            fwi = atof(token);  token=strtok(NULL, ",");
            strcpy(fbp.FuelType, token);  token=strtok(NULL, ",");
                if (!strcasecmp(fbp.FuelType, "O1") ) strcpy(fbp.FuelType, "O1b");
            
            ros = atof(token);  token=strtok(NULL, ",");
            sfc = atof(token);  token=strtok(NULL, ",");
            tfc = atof(token);  token=strtok(NULL, ",");
            bfc = atof(token);  token=strtok(NULL, ",");
            hfi = atof(token);  token=strtok(NULL, ",");
            
// ... estarea,UTC,...
             estarea = atof(token);  token=strtok(NULL, ",");
            strcpy(cmc_date, token);  token=strtok(NULL, ",");

if (feps.print<0 && strlen(cmc_date)<11) printf("2.b. Warning: incomplete cmc_date (%s),\n", cmc_date);
            
/*  # temp = Surface temperature (degree Celcius)
    # rh = Surface relative humidity (fraction - unitless)
    # ZS =  Surface geopotential height (dam)
    # ws-met =  Surface wind speed (knots)
    # Precip rate = Total precipitation rate (m/s)
    # TS = Surface and soil temperatures (Superficial) (degree Kelvin)
    # T850 = Temperature @ 850 mb (degree Celcius)
    # T700 = Temperature @ 700 mb (degree Celcius)
    # T500 = Temperature @ 500 mb (degree Celcius)
    # Z850 = Geopotential height @ 850 mb (dam)
    # Z700 = Geopotential height @ 700 mb (dam)
    # Z500 = Geopotential height @ 500 mb (dam) */

            feps.temp = atof(token);  token=strtok(NULL, ",");
            feps.rh = atof(token);  token=strtok(NULL, ",");
            feps.zs = atof(token);  token=strtok(NULL, ",");
            fbp.WS = atof(token);  token=strtok(NULL, ",");
            precip = atof(token);  token=strtok(NULL, ",");
            feps.ts = atof(token);  token=strtok(NULL, ",");
            feps.t850 = atof(token);  token=strtok(NULL, ",");
            feps.t700 = atof(token);  token=strtok(NULL, ",");
            feps.t500 = atof(token);  token=strtok(NULL, ",");
            feps.t250 = atof(token);  token=strtok(NULL, ",");
            feps.z850 = atof(token);  token=strtok(NULL, ",");
            feps.z700 = atof(token);  token=strtok(NULL, ",");
            feps.z500 = atof(token);  token=strtok(NULL, ",");
            feps.z250 = atof(token);

// some values may be reported in Kelvin
            if (feps.temp > 100) feps.temp = feps.temp - 273.16;
            if (feps.ts > 100) feps.ts = feps.ts - 273.16;
            if (feps.t850 > 100) feps.t850 = feps.t850 - 273.16;
            if (feps.t700 > 100) feps.t700 = feps.t700 - 273.16;
            if (feps.t500 > 100) feps.t500 = feps.t500 - 273.16;
            if (feps.t250 > 100) feps.t250 = feps.t250 - 273.16;

// in case RH is reported as a fraction
            if (feps.rh < 1.) feps.rh = 100.*feps.rh;
            
            fbp.WS = 1.852 * fbp.WS;  /* knots to km/hr */
            precip = 3600000.* precip;  /* m/s to mm/hr */

// in case heights are reported in decametres
            if (feps.z850 < 500) feps.zs = 10* feps.zs;  // in the case of the surface, assume that the units are tyhe same as those for 850 mb
            if (feps.z850 < 500) feps.z850 = 10* feps.z850;
            if (feps.z700 < 1000) feps.z700 = 10* feps.z700;
            if (feps.z500 < 1000) feps.z500 = 10* feps.z500;
            if (feps.z250 < 2000) feps.z250 = 10* feps.z250;
            
            feps.zt = feps.z250;    // arbitrarily setting the plume top to 250 mb

/* 2.c. Assign upper air values into arrays uP, uZ and uT for future use in PlumeCalcDry (under development) */
// if (feps->print==-2)printf("2.c.\n");

            ua.iLevels=0;
            
            ua.P[0] = 1000.;     // should use surface pressure but am using 1000 mb for now
            
            ua.P[0] = 850.-(feps.z850-feps.zs)/(feps.z700-feps.z850)*(700.-850.); 
            ua.Z[0] = feps.zs;
            ua.T[0] = feps.ts;

            ua.iLevels++;

            if (feps.z850 > feps.zs)
            {
                ua.P[ua.iLevels] = 850.;
                ua.Z[ua.iLevels] = feps.z850;
                ua.T[ua.iLevels] = feps.t850;
                ua.iLevels++;
            }

            if (feps.z700 > feps.zs)
            {
                ua.P[ua.iLevels] = 700.;
                ua.Z[ua.iLevels] = feps.z700;
                ua.T[ua.iLevels] = feps.t700;
                ua.iLevels++;
            }

            ua.P[ua.iLevels] = 500.;
            ua.Z[ua.iLevels] = feps.z500;
            ua.T[ua.iLevels] = feps.t500;
            ua.iLevels++;

            ua.P[ua.iLevels] = 250.;
            ua.Z[ua.iLevels] = feps.z250;
            ua.T[ua.iLevels] = feps.t250;
            ua.iLevels++;

            ua.P[ua.iLevels] = 100.;
            ua.Z[ua.iLevels] = feps.z250*2;
            ua.T[ua.iLevels] = feps.t250;    // assuming isothermic above 250 mb
            ua.iLevels++;
 
/* 2.d. Determine if the line read in is a new hotspot 
        (if the latitude, longitude or date does not equal previous line’s values).  
        If new, reset arrays to zero. */ 
// if (feps->print==-2)printf("2.d.\n");

            if (plat==fbp.LAT && plon==fbp.LON && cmcdate <= atoi(&cmc_date[0])) // same fire
            {
                newfire = 0; // false
                // feps.area=parea;        // overwrite area with area from previous time step
            }
            else    // new fire
            {
                newfire = 1; // true
                feps.area=0.; //=estarea;
                totalemissions = 0.;  // zero everything including total emissions
                feps.Qo = Qo = 0.;    // Qo is used to temporarily store the energy when comparing height predictions using different lapse rates 
                feps.Qplume = 0;
                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
                    for (iUTC=0; iUTC<nlong+1; iUTC++)
                    {
                        Qs[iUTC]=0.;
                        emissions.f[iSpecies][iUTC]=0.;       // new fire, set concetrations back to zero
                        emissions.s[iSpecies][iUTC]=0.;       // new fire, set concetrations back to zero
                        emissions.r[iSpecies][iUTC]=0.;       // new fire, set concetrations back to zero
                    }
                plat=fbp.LAT;
                plon=fbp.LON;
                parea=0.;

                feps.growth=0.;
                feps.Flaming=0.;
                feps.Smoldering=0.;
                feps.Residual=0.;
                feps.dz = dz = -999.;
                dznew = -999.;
                dznew = 0.;
                
                for (i = 0; i < MAX_HOURS+1; i++)
                {
                    Fs[i]=0.;       // new fire, set concentrations back to zero
                    Ss[i]=0.;       // new fire, set concentrations back to zero
                    Rs[i]=0.;       // new fire, set concentrations back to zero
                }
            }

/* 2.e. Reassign FBP fuel type (if required).  If there is no matching FBP fuel type, set to NF. */            
// if (feps->print==-2)printf("2.e.\n");
            if (!strcasecmp(fbp.FuelType,"cropland")) strcpy(fbp.FuelType, "O1a");
            if (!strcasecmp(fbp.FuelType,"low_veg")) strcpy(fbp.FuelType, "O1a");  
            if (!strcasecmp(fbp.FuelType,"urban")) strcpy(fbp.FuelType, "NF");
            if (!strcasecmp(fbp.FuelType,"bog")) strcpy(fbp.FuelType, "NF");
            if (!strcasecmp(fbp.FuelType,"water")) strcpy(fbp.FuelType, "NF");
            if (!strcasecmp(fbp.FuelType,"non_fuel")) strcpy(fbp.FuelType, "NF");
            if (!strncasecmp(fbp.FuelType,"M1", 2)) strcpy(fbp.FuelType, "M1"); 
            if (!strcasecmp(fbp.FuelType,"D2")) strcpy(fbp.FuelType, "D1");
            
/* if there is no matching FuelType to FBP, set to NF */
            if (! (    !strncasecmp(fbp.FuelType, "C1", 2) ||
                    !strncasecmp(fbp.FuelType, "C2", 2) ||
                    !strncasecmp(fbp.FuelType, "C3", 2) ||
                    !strncasecmp(fbp.FuelType, "C4", 2) ||
                    !strncasecmp(fbp.FuelType, "C5", 2) ||
                    !strncasecmp(fbp.FuelType, "C6", 2) ||
                    !strncasecmp(fbp.FuelType, "C7", 2) ||
                    !strncasecmp(fbp.FuelType, "D1", 2) ||
                    !strncasecmp(fbp.FuelType, "D2", 2) ||
                    !strncasecmp(fbp.FuelType, "M1", 2) ||
                    !strncasecmp(fbp.FuelType, "M2", 2) ||
                    !strncasecmp(fbp.FuelType, "M3", 2) ||
                    !strncasecmp(fbp.FuelType, "M4", 2) ||
                    !strncasecmp(fbp.FuelType, "S1", 2) ||
                    !strncasecmp(fbp.FuelType, "S2", 2) ||
                    !strncasecmp(fbp.FuelType, "S3", 2) ||
                    !strncasecmp(fbp.FuelType, "O1", 2) ) )
                strcpy(fbp.FuelType, "NF");
            
/* 2.f. Determine time zone (longitude/15) and parse detection date and time.  
        Adjust detection time from UTC to local time such that noon is at solar zenith . 
        Adjust LDT to LST (hardcoded to LST?).  
        Determine UTC and LST time of the forecast record. */            
// if (feps->print==-2)printf("2.f.\n");
            timezone = (fbp.LON)/15.;  
            feps.LDT = 0.;

// e.g. 20160216_00:00
            cmcdate = atoi(&cmc_date[0]);
            cmcUTC = atoi(&cmc_date[9]);
            cmcLST = cmcUTC + (int)timezone;
            if (cmcLST < 0) cmcLST = cmcLST+24;

/* Determine CMC julian date */
            if (cmcdate> 0)
            {
                
                if (cmcdate > 9999)
                    iyear = (int)(cmcdate/10000);
                else
                    iyear = 0;
                imonth = (int)((cmcdate-10000*iyear)/100);
                iday = cmcdate-10000*iyear-100*imonth;
                cmcDj = jdays[imonth-1]+iday;
                if (iyear > 0 && iyear%4 == 0 && imonth >2) cmcDj+=1;

            }  // iyear, imonth and iday will be overwritten

/* Determine detection date and time */
            // 2014-04-01 20:01:00
            feps.ddate = atoi(&rep_date[0])*10000 + atoi(&rep_date[5])*100 + atoi(&rep_date[8]);
            feps.dtime  = atoi(&rep_date[11])*100 + atoi(&rep_date[14]);    // currently in UTC
            feps.timestep = 1.;

/* Determine julian date */
            if (fbp.Dj = -1)
            {
                if (feps.ddate > 9999)
                    iyear = (int)(feps.ddate/10000);
                else
                    iyear = 0;
                imonth = (int)((feps.ddate-10000*iyear)/100);
                iday = feps.ddate-10000*iyear-100*imonth;
                fbp.Dj = jdays[imonth-1]+iday;
                if (iyear > 0 && iyear%4 == 0 && imonth >2) fbp.Dj+=1;
            }

/* Adjust LDT value to LST */
            if (feps.LDT == 100 || feps.LDT == 1)
                feps.LDT = 1;
            else
                feps.LDT = 0;

            detectionUTC = (int)(feps.dtime / 100) + (feps.dtime % 100) /60.;
            detectionLST = (int)(feps.dtime / 100) + (feps.dtime % 100) /60. - feps.LDT + (int)timezone;

            if (detectionLST < 0.) 
            {
                detectionLST = detectionLST + 24.;
                fbp.Dj = fbp.Dj;    // keep fbp.Dj as the UTC date
            }


/* 2.g. Recalculate the BUI value from DMC and DC values.  Determine whether the BUI effect is being used. */    
// if (feps->print==-2)printf("2.g.\n");
            fbp.BUI = BUIcalc(feps.DMC, feps.DC);

/* If BUI is zero or less, BUIeff = False */
            if (fbp.BUI <= 0) fbp.BUIEff = -1;

/* 2.h. If new fire, determine the hours from the first hourly CMC forecast to the detection time: 
           if less than 24 hours, an initial fire size is calculated (based on persistence routine); 
           if greater than zero, then the fire starts at zero size. */
// if (feps->print==-2)printf("2.h.\n");

            if (newfire == 1) // determine the initial Area(t) for first CMC hour.
            {
                dt = 24.*(fbp.Dj-cmcDj) + detectionUTC + cmcUTC; // hours from first CMC hour to detection time
                
                if (dt < 24.)
                {
                    
                    if (!strncmp(feps.shape, "tophat", 3) ||
                        !strncmp(feps.shape, "weighted", 3))
                    {
                        
/* new approach, count backwards from detectionLST */
                        sum = weight[(int)detectionLST] * (detectionLST - (int)detectionLST);
                
                        if (cmcLST>detectionLST) // cmcLST on previous day
                            for (i=(int)detectionLST+24; i>cmcLST;i--)
                                if (i>23)
                                    sum+=weight[i-24]; // accumulate the hourly weight over the course of a day
                                else
                                    sum+=weight[i]; // accumulate the hourly weight over the course of a day
                        else
                        {
                            for (i=(int)detectionLST; i>cmcLST;i--)
                            {
                                if (i>23)
                                    sum+=weight[i-24]; // accumulate the hourly weight over the course of a day
                                else
                                    sum+=weight[i]; // accumulate the hourly weight over the course of a day
                            }
                        }

                        if (dt<0)
                            estarea2 = estarea - estarea*sum + (1-(int)dt/24)*estarea;  // account for burning prior to forecast start
                        else                
                            estarea2 = estarea - estarea*sum;  // new value for estarea, capturing the entire burn period

                    }

                }
                else        // detection 24 hours or more after first CMC hour
                    estarea2 = 0.;
                    
                feps.area = parea = estarea2;
            }
            


/* adjust cmcUTC to a time wrt detection time and date, so this is 0 at detection time, negative before and positive after (can exceed 24) */
            cmcUTC = cmcUTC + 24*(cmcDj - fbp.Dj);    // should use a new variable name like "elapsed"
            if (newfire == 1) cmcUTCo = cmcUTC; // true, cmcUTCo = the first fire record time (not the time of detection)

if (feps.print<0 && cmcUTC<cmcUTCo) printf("2.h. Warning: cmcUTC (%d) < cmcUTCo (%d)", cmcUTC, cmcUTCo);

/* elapsed is the time since detectionLST (dhrs) */
            elapsed = 0.;firsttime = 1;
            timestep = feps.timestep;

/* 2.i. Calculate FFMC from yesterday’s value using the technique described in Lawson et al. (1996).  
        This is based on the observed noon FFMC.  
        If the time is before noon, yesterday’s FFMC is required.  
        This is converged on so that yesterday’s FFMC will predict today’s observed FFMC 
        (set to 11:59AM, 23:59 hours later). */
// if (feps->print==-2)printf("2.i.\n");

            FFMCy = FFMC = fbp.FFMC;     /* we are assuming the value passed in FBP is the noon FFMC */

 /* Calculate yesterday's FFMC */
            if (detectionLST < 12.)
            {
                FFMCy = FFMC; FFMCt=0.;    /* first guess */
                FFMCnoon = bdl_ffmc(FFMC, 1200, (int)feps.rh);
                i=0;
                while (fabs(FFMCnoon-FFMCt) > .1 && FFMCy > 50. && FFMCy < 100. && i < 100)
                {
                    FFMCt = bdl_ffmc(FFMCy, 1159, (int)feps.rh);
                    FFMCy = FFMCy + (FFMCnoon-FFMCt);
                    i++;
                }
                if (i==100) FFMCy = (FFMCt+FFMCy)/2;
                if (FFMCy>99.) FFMCy=99.;
            }

            if (FFMCy < 60. && strcasecmp(feps.method, "AB") && strcasecmp(feps.method, "cmc") )    // we're already in !strcasecmp(feps.method, "cmc") 
            {
                fprintf(outfile, "Low FFMC: %lf set to 75.\n", FFMCy);
                FFMCy = 75.;
            }

            hour = mtime(detectionLST);

            if (hour < 1200)
                fbp.FFMC = bdl_ffmc(FFMCy, hour, (int)feps.rh);
            else
                fbp.FFMC = bdl_ffmc(FFMC, hour, (int)feps.rh);
            fbp.t = detectionLST;

/* 2.j  Calculate FBP values at time of detection (and then at current time) */
// if (feps->print==-2)printf("2.j.\n");
            status = FBPCalc (&fbp);    // FBP at detectionLST

/* 2.k. Determine discovery time (the number of hours of growth required to reach the detection size).  
        When persistence is being used, this time is assumed to be 24 hours. */
// if (feps->print==-2)printf("2.k.\n");
            t1=0.; t2=0.; A=estarea; P=0.; D=0.;
            PROCalc(fbp.FuelType, fbp.Accel, 
                        fbp.CFB, fbp.ROS, fbp.FROS, fbp.BROS,   // ROS at UTC time
                        &t1, &t2, &A, &P, &D);

            discovery = t2;    // use the detectionLST to determine discovery; afterwards use cmcLST for the hour 

            if (!strncmp(feps.shape, "tophat", 3) ||
                !strncmp(feps.shape, "weighted", 3)) 
                    discovery = 24.; // set the time of discovery equal to 24 hours 

            hour = (int)(100*cmcLST);
            feps.elapsed = cmcUTC - detectionUTC;  // used in Brigg's model to determine the current time
            if (feps.elapsed < 0.) feps.elapsed = 0.;
            
            if (hour < 1200)
                fbp.FFMC = bdl_ffmc(FFMCy, hour, (int)feps.rh);
            else
                fbp.FFMC = bdl_ffmc(FFMC, hour, (int)feps.rh);

            
/* 2.j  Calculate FBP values at current time */
            status = FBPCalc (&fbp);    // FBP at UTC time

            if (status==0)
            {
                
/* 2.l. If current time is prior to discovery, the fire has not yet occurred and does not grow; else: */
// if (feps->print==-2)printf("2.l.\n");
//                if ( ( (detectionUTC-discovery) >= (double)cmcUTC ) || discovery == 999.9 ||  FFMCy < 60.)  // no fire yet
                if ( ( (detectionUTC-discovery) >= (double)cmcUTC ) || discovery == 999.9)  // no fire yet
                {
                    A = 0.;
                    feps.area = 0.;
                    r_smoke = -999.;
                    r_smoke = 0.;
                    dz = -999;
                    dH = -999.;
                    discovery = 0;
                }
                else 
                {

/* 2.l.i. If the current time corresponds to the reset time (if being used, not the default), set the  area, cumulative energy and total emissions to zero */
// if (feps->print==-2)printf("2.l.i\n");
                    if ( (feps.reset >= 0.) && (cmcLST >= feps.reset && cmcLST < feps.reset+timestep) )
                    {
                        parea = 0.;
                        Qo=0.;
                        totalemissions = 0.;
                    } 
                    else
                    {
                        parea = feps.area;
                        Qo=feps.Qo;
                    }
/* 2.l.ii. Determine difference in time between current time and detection time.  If less than an hour, calculate partial hour growth; otherwise, calculate hourly growth according to top-hat or weighted approach.  Add to total area. */            
// if (feps->print==-2)printf("2.l.ii.\n");
                    if (!strncmp(feps.shape, "elliptical", 3))    // this portion of the code seems to be in flux
                        if (parea ==  0.)
                        {
                            t1 = 0;
                            t2 = 24.*( (fbp.Dj+(double)cmcUTC/24.) - (fbp.Dj+(detectionUTC-discovery)/24.) );    // these can be multiday lengths, when grown at noon ROS produce wrong values
                            t2 = (double)cmcUTC - (detectionUTC-discovery) ;    // these can be multiday lengths, when grown at noon ROS produce wrong values

                            if (t2 > 0.)
                            {
                                A=0.; P=0.; D=0;
                                PROCalc(fbp.FuelType, fbp.Accel, 
                                        fbp.CFB, fbp.ROS, fbp.FROS, fbp.BROS,
                                        &t1, &t2, &A, &P, &D);
                                feps.growth = A;   // area at cmcUTC
                            }
                            else
                                feps.growth = 0.;

                            if (cmcUTCo >=0) parea = estarea;

                            feps.growth = estarea;  //temporary fix
                        
                            parea = estarea2; // new approach 2018/01/18
                            feps.growth = 0;  //temporary fix

                        }   // if (parea ==  0.)                   
                        else
                        {
                            t1=0.; t2=0.; A=parea; P=0.; D=0.;                         
                            PROCalc(fbp.FuelType, fbp.Accel, 
                                    fbp.CFB, fbp.ROS, fbp.FROS, fbp.BROS,
                                    &t1, &t2, &A, &P, &D);

                            t1=t2; t2=t2+timestep; A=0.;P=0.;D=0.;
                            PROCalc(fbp.FuelType, fbp.Accel, 
                                    fbp.CFB, fbp.ROS, fbp.FROS, fbp.BROS,
                                    &t1, &t2, &A, &P, &D);

                            feps.growth = A;  // area growth in one time step: A(t2+timestep)-A(t2)
                        }

                    dt = (double)cmcUTC - (detectionUTC-discovery); // difference in time between cmcUTC and detectionUTC

                    if (!strncmp(feps.shape, "tophat", 3))
                    {
                        if (hour > mtime(feps.thstart) && hour <= mtime(feps.thend))
                            if (dt > 0.  &&  dt < 1.)
                                feps.growth = dt * timestep * estarea/(feps.thend - feps.thstart);      // top-hat
                            else
                                feps.growth = timestep * estarea/(feps.thend - feps.thstart);      // top-hat
                        else
                            feps.growth = 0.;

                        
/*                        if (hour > 900 && hour <= 2100)
                            if (dt > 0.  &&  dt < 1.)
                                feps.growth = dt * timestep * estarea/12.;      // top-hat
                            else
                                feps.growth = timestep * estarea/12.;      // top-hat
                        else
                            feps.growth = 0.; */
                    } 
                    
                    if (!strncmp(feps.shape, "weighted", 3))    // this is the default based on readfeps()
                    {
                        if (cmcLST >= 24) cmcLST-=24;
                        
                        if (dt > 0.  &&  dt < 1.)
                            feps.growth = dt* weight[cmcLST] * estarea;      // weighted fire growth of previous hour 
                        else
                            feps.growth = weight[cmcLST] * estarea;      // weighted fire growth of previous hour 

                    }
                    
                    if (!strncmp(feps.shape, "line", 3) )
                    {   
                        if (parea == 0. && dt > 0.) parea = estarea;
                            
                        diameter = 2*sqrt(10000.*estarea/pi);    // diameter of a circle of estarea

                        if (dt > 0.  &&  dt < 1.)
                            feps.growth = timestep * dt * 60.*fbp.ROS * diameter/10000.;
                        else
                            feps.growth = timestep * 60.*fbp.ROS * diameter/10000.;
                    }
                    
                    if (newfire==1) feps.growth = 0.;
                    
                    feps.area = feps.growth + parea;
                    parea= feps.area;

                    dz8 = dz7 = dz5 = dz = -999.;
                    l8 = l7 = l5 = -999.;
                    r_smoke = -999.;
                    r_smoke = 0.;

/* 2.l.iii. Calculate energy of the fire */
// if (feps->print==-2)printf("2.l.iii.\n");
                    feps.Qplume = EnergyCalc(&feps, &fbp);

/* 2.l.iv. Calculate the energy released from flaming, smoldering and residual combustion.
            Calculate flaming, smoldering and residual emissions rates in tonnes/h  */                    
// if (feps->print==-2)printf("2.l.iv.\n");
                    EmissionsCalc(&feps, &fbp);   
                        emissions.residence_flaming=0.25;   /* residence time for flaming combustion (hours) -- remains the same */
                        emissions.residence_smoldering = feps.DOB_S/2.; /* this is new, based on the assumption that a fire smolders at 1 cm/hour (Dan Thompson) */
                        emissions.residence_residual = feps.DOB_S/2.;   /* split equally between smoldering and residual */

                    jModel = -1;    /* an admittedly awkward way of creating hourly energy components */

/* 2.l.v. Distribute the energy over time */                    
// if (feps->print==-2)printf("2.l.v.\n");
                    iMax = EmissionsOverTime(Fs, Ss, Rs, jModel, iSpecies, &emissions, &feps);

                    for (i=0; i< iMax; i++)
                        if (i+cmcUTC-cmcUTCo <= MAX_HOURS)
                        {
                            Qs[i+cmcUTC-cmcUTCo] = Qs[i+cmcUTC-cmcUTCo] + Fs[i]+Ss[i]+Rs[i];
                        }
                        else 
                            if (feps.print==-1) printf("2.l.v. Warning: i+cmcUTC-cmcUTCo (%d) . MAX_HOURS (%d)\n", i+cmcUTC-cmcUTCo, MAX_HOURS);

/* 2.l.vi. From the area growth and energy released, calculate the plume rise */                        
// if (feps->print==-2)printf("2.l.vi.\n");
                    feps.Qplume = Qs[cmcUTC-cmcUTCo];  // need to check this 
// printf("%s %e %e %e\n", cmc_date, feps.Qplume, feps.Qo, Qo);
                    feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                
                    dznew = PlumeCalcDry(&feps, &ua); // new piece-wise approach -- included here for comparative purposes to the standard

/* 2.l.vi. (1) if type=dry, use all upper air data following a piecewise method of calculating plume rise */                
                    if (!strcmp(feps.type, "dry"))
                    {
                        feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                        dz = PlumeCalcDry(&feps, &ua); // new piece-wise approach
                    }
/* 2.l.vi. (2) follow the standard average lapse rate method outlined in Anderson et al., 2011 */
                    else    // follow the standard average lapse rate method outlined in Anderson et al., 2011
                    {

/* 2.l.vi. (2) (a) calculate height based on surface to 850 mb lapse */
                        feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                        l8 = feps.lapse = (feps.t850 - feps.ts) / (feps.z850 - feps.zs);
                        if (feps.lapse > -0.0098 && feps.z850 > feps.zs)
                            dz8 = dz = PlumeCalc(&feps);

/* 2.l.vi. (2) (b) if recalculated height greater than 2000 m, 
                  recalculate height based on surface to 700 mb lapse rate
                  (if resulting height less than 2000, use average) */
                        if (dz > 2000.  || dz <= 0.)
                        {
                            feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                            l7 = feps.lapse = (feps.t700 - feps.ts) / (feps.z700 - feps.zs);
                            if (feps.lapse > -0.0098 && feps.z700 > feps.zs)
                            {
                                dz7 = dz = PlumeCalc(&feps);
                                if (dz < 2000)
                                {
                                    feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                                    feps.lapse = (l8+l7)/2;
                                    if (feps.lapse > -0.0098)
                                        dz = PlumeCalc(&feps);
                                }
                            }
                
/* 2.l.vi. (2) (c) if recalculated height greater than 4000 m, 
                  recalculate height based on surface to 500 mb lapse rate  
                  (if resulting height less than 4000, use average) */
                            if (dz > 4000.  || dz <= 0.)
                            {
                                feps.Qo = Qo;  // remember the energy prior to entering decision tree 
                                l5 = feps.lapse = (feps.t500 - feps.ts) / (feps.z500 - feps.zs);
                                if (feps.lapse > -0.0098)
                                {
                                    dz5 = dz = PlumeCalc(&feps);
                                    if (dz < 2000)
                                    {
                                        feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                                        feps.lapse = (l8+l7)/2;
                                        if (feps.lapse > -0.0098)
                                            dz = PlumeCalc(&feps);                        
                                    }       
                                }
                                
/* 2.l.vi. (2) (d) if superadiabatic from surface up to 500m, recalculate using 850 mb as the base
                  recalculate height based on 850 mb to 700 mb lapse rate
                  (if resulting height less than 2000, use average) */
                                if (dz <= 0.)
                                {
                                    feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                                    l7 = feps.lapse = (feps.t700 - feps.t850) / (feps.z700 - feps.z850);
                                    if (feps.lapse > -0.0098  && feps.z700 > feps.zs)
                                    {
                                        dz7 = dz = PlumeCalc(&feps);
                                        if (dz < 2000)
                                        {
                                            feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                                            feps.lapse = (l8+l7)/2;
                                            if (feps.lapse > -0.0098)
                                                dz = PlumeCalc(&feps);
                                        }
                                    }
                
/* 2.l.vi. (2) (e) if recalculated height greater than 4000 m, 
                  recalculate height based on 850 mb to 500 mb lapse rate  
                  (if resulting height less than 4000, use average) */
                                    if (dz > 4000.  || dz <= 0.)
                                    {
                                        feps.Qo = Qo;  // remember the energy prior to entering decision tree 
                                        l5 = feps.lapse = (feps.t500 - feps.t850) / (feps.z500 - feps.z850);
                                        if (feps.lapse > -0.0098)
                                        {
                                            dz5 = dz = PlumeCalc(&feps);
                                            if (dz < 2000)
                                            {
                                                feps.Qo = Qo;   // remember the energy prior to entering decision tree 
                                                feps.lapse = (l8+l7)/2;
                                                if (feps.lapse > -0.0098)
                                                    dz = PlumeCalc(&feps);                        
                                            }       
                                        }
                                    }    // if (dz > 4000.  || dz < 0.) 
                                }    // if (dz < 0.)  [superadiabatic from surface]        
                            }    // if (dz > 4000.  || dz < 0.) 
                        }    // if (dz > 2000.  || dz < 0.)           
                    }    // else 

                    if (dz8 > 2000 && dz7 < 2000) 
                    {
                        dz = (dz8+dz7)/2;
                        feps.lapse = (l8+l7)/2;
                    }

                    if (dz7 > 4000 && dz5 < 4000) 
                    {
                        dz = (dz7+dz5)/2;
                        feps.lapse = (l7+l5)/2;
                    }

/* 2.l.vi. (3) limit height to 250 mb height */
                    if (dz > feps.z250)    dz = feps.z250;

/* 2.l.vi. (4) impose a minimum height */
                    if (dz < feps.dzmin) dz = feps.dzmin;

                    feps.Qo = Qo;  // remember the energy prior to entering decision tree (note that as is, this discounts the energy lost due to black body radiation)

                    feps.Qo = Qo = Qo + feps.Qplume;    // I used to do this in PlumeCalc but the new EmissionsCalc has thrown this off

                    feps.dz = dz; 

/* The Briggs calculations are highly unreliable */                    
//                    dH = OldBriggsCalc(&feps);
//                    dH = BriggsCalc(&feps);

                } // if ( (fbp.Dj+(detectionUTC-discovery)/24.) > (cmcDj+(double)cmcUTC/24.)  || (discovery == 999.9) ||  FFMCy < 60.)  ...  else 

/* 2.l.vii. Find emodel that matches the fuel type */                
// if (feps->print==-2)printf("2.l.vii.\n");
                if (emissions.nModels==0 || emodels.nModels==0)
                    jModel=0;
                else
                    for (iFuel=0; iFuel<nFuels; iFuel++)
                    if (!strcasecmp(emodels.FuelType[iFuel], fbp.FuelType)) // found a match
                        for (iModel=0;iModel<emissions.nModels;iModel++)
                            if (emissions.emodel[iModel]==emodels.emodel[iFuel] )
                                jModel=iModel; 

/* 2.l.viii.  Distribute the emission over time */
// if (feps->print==-2)printf("2.l.viii.\n");
                for (iSpecies = 0; iSpecies < emissions.nSpecies; iSpecies++)
                {   
                    iMax = EmissionsOverTime(Fs, Ss, Rs, jModel, iSpecies, &emissions, &feps);
                    for (i=0; i< iMax; i++)
                        if (i+cmcUTC-cmcUTCo <= MAX_HOURS)
                        {
                            emissions.f[iSpecies][i+cmcUTC-cmcUTCo]=emissions.f[iSpecies][i+cmcUTC-cmcUTCo]+Fs[i];
                            emissions.s[iSpecies][i+cmcUTC-cmcUTCo]=emissions.s[iSpecies][i+cmcUTC-cmcUTCo]+Ss[i];
                            emissions.r[iSpecies][i+cmcUTC-cmcUTCo]=emissions.r[iSpecies][i+cmcUTC-cmcUTCo]+Rs[i];
                        }
                        else
                            if (feps.print==-1) printf("2.l.viii. Warning: i+cmcUTC-cmcUTCo (%d) . MAX_HOURS (%d)\n", i+cmcUTC-cmcUTCo, MAX_HOURS);

                }

/* 2.l.ix. Calculate total emissions and r_smoke */
// if (feps->print==-2)printf("2.l.ix.\n");

                totalemissions = totalemissions+10.*fbp.TFC*feps.growth;

                if (feps.dz > 0 && feps.M > 0) 
                    r_smoke = 1000 * (1000.* totalemissions) / feps.M;
                else
                    r_smoke = 0.;
                

/* replace missing values (-9999) or any negative plume heights with 0. */
                if (feps.dz < 0.)
                {
                    feps.dz = 0.;
                    feps.M = 0.;
                }


            }  // if (strcasecmp(fbp.FuelType, "NF"))
            else
            {
                A = 0.;
                feps.area = 0.;
                feps.growth = 0.;
                totalemissions = 0.;
                r_smoke = 0.;
                feps.dz = 0;
                feps.M = 0.;
                dH = -999.;
                discovery = 0;
                feps.Qo = 0.;             
                feps.Qplume = 0.; 
                dznew = 0.;
            } 
            
/* 2.m.  Print results to the screen. */            
// if (feps->print==-2)printf("2.m.\n");
            
// if (cmcUTC==24) // a trick that prints out just the 24 hour growth (Area(t) - estarea), useful for evaulation
            {
// lat, lon, rep_date, source, sensor, ffmc, dmc, dc, ws, fwi, fuel, ros, sfc, tfc, bfc, hfi, estarea, 
                fprintf(outfile, "%.3lf, %.2lf, %s, %s, %s, %.1lf, %.1lf, %.1lf, %.1lf, %.1lf, %s, %lf, %lf, %lf, %lf, %.0lf, %.2lf,", 
                    fbp.LAT, fbp.LON, rep_date, source, sensor, fbp.FFMC, feps.DMC, feps.DC, fbp.WS, 
                    fwi, fbp.FuelType, fbp.ROS, fbp.SFC, fbp.TFC,bfc, fbp.HFI, estarea);

                fprintf(outfile, " %s, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf, %.2lf,", 
                    cmc_date, feps.temp, feps.rh, fbp.WS, precip, feps.ts, feps.t850, feps.t700, feps.t500, feps.t250, feps.zs, feps.z850, feps.z700, feps.z500, feps.z250);

//                fprintf(outfile, " %lf, %lf, %lf, %lf, %d, %d", feps.area, feps.growth, totalemissions, r_smoke, (int)feps.dz, (int)dznew); 
                fprintf(outfile, " %lf, %lf, %lf, %lf, %d, %lf", feps.area, feps.growth, totalemissions, r_smoke, (int)feps.dz, feps.M); 
                fprintf(outfile, ", %lf, %lf", feps.Qo, feps.Qplume);

                for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
                {
                    fprintf(outfile, ", %.4e", emissions.f[iSpecies][cmcUTC-cmcUTCo]);
                    fprintf(outfile, ", %.4e", emissions.s[iSpecies][cmcUTC-cmcUTCo]);
                    fprintf(outfile, ", %.4e", emissions.r[iSpecies][cmcUTC-cmcUTCo]);
                }
                fprintf(outfile, "\n");
            } // just track 24 hour area and growth (without the "if (cmcUTC==24)" statement, the brackets have no effect)
     

/* 2.n. Read next line */
// if (feps->print==-2)printf("2.n.\n");
            fgets(buffer, 480, csvfile);

        } while (!feof(csvfile) );  

    }  //   if (!strcasecmp(feps.method, "cmc") )  /* CMC Firework method *////////////////////////////////////////////////////////////////////////////////////////////
    
/* the following section is included to test PlumeCalcDry(), the piece-wise plume rise calculation, following an ICAO standard atmosphere */

    else if (!strcasecmp(feps.method, "ICAO") )
    {
        feps.Qo = 0.;
        feps.Qplume = 8.35E10*feps.area;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m AGL

// from data stored in infile.dat
        feps.dz = PlumeCalc(&feps);
        printf("PlumeCalc (infile.dat): %lf\n\n", feps.dz);    

//         printf("Hello: %s\n\n", feps.method);

        /* ICAO Standard Atmosphere */

        feps.ts = 14.3;
        feps.zs = 111.;
        feps.pressure = 1000.;
        feps.area = 1.0;
        feps.perimeter=0.;
        feps.residencetime=0.;
        feps.radiation = 100;
        feps.alpha=0.;
        feps.lapse = -6.5/1000.;
        feps.growth=-1.;
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m AGL
        feps.dz = PlumeCalc(&feps);
        printf("PlumeCalc (lapse rate): %lf\n\n", feps.dz);

        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.lapse = (feps.t500-feps.t850)/(feps.z500-feps.z850);
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalc(&feps);
        printf("PlumeCalc (%s): %lf\n\n", profilename, feps.dz);        

        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.lapse = -999.;
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalc(&feps);
        printf("PlumeCalc (%s): %lf\n\n", profilename, feps.dz);        

printf("**************************\n");        
        
        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalcDry(&feps, &ua);    
        printf("PlumeCalcDry (%s): %lf\n\n", profilename, feps.dz);        

        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm_8752.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalcDry(&feps, &ua);    
        printf("PlumeCalcDry (%s): %lf\n\n", profilename, feps.dz);        

        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm_752.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalcDry(&feps, &ua);    
        printf("PlumeCalcDry (%s): %lf\n\n", profilename, feps.dz);        

        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm_852.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.Qo = 0.;
        feps.Qplume = 1636.5*5102.*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalcDry(&feps, &ua);    
        printf("PlumeCalcDry (%s): %lf\n\n", profilename, feps.dz);        

        sprintf(profilename, "D:\\Alberta Smoke Monitoring\\Soundings\\ICAO_2015062718_gpm_8752_var.profile");
        iLevels = Readupperair(profilename, &feps, &ua);
        feps.Qo = 0.;
        feps.Qplume = (1636.5+254.8364)*5102*10000.*feps.area;    // using a standard atmosphere, plume height should be 5463 m
        feps.dz = PlumeCalcDry(&feps, &ua);    
        printf("PlumeCalcDry (%s): %lf\n\n", profilename, feps.dz);        

// feps.Qplume = 1842.881*5102*10000.*feps.area;   // using ICAO2 with an isothermal layer under 800 mb, should get the same height using this energy
// feps.Qplume = (1842.881+254.8364)*5102*10000.*feps.area;
// feps.Qplume = 2097.717*5102*10000.*feps.area;        
    }
        
    else    // the old code (removed)

    {  

    }// ! //   if (!strcasecmp(feps.method, "cmc") )  ... else /* CMC Firework method */

    return (status);
}

