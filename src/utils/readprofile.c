#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include "thermo.h"
#include "fbp_2009.h"
#include "feps_2011.h"
//#include "malloc.h"

double rh_calc();

FILE *infile;


double Readprofile(char *infilename, FEPS *feps, FBP *fbp)
{
char header[300], buffer[300], id[25], *args[10], *token;
double lat, lon, start_time, PREV_PRES, PREV_TEMP,
        PRSS, TPP1, SHGT, MXHT, SHTF, USTR, DSWF,
        UWND, VWND, TEMP, SPHU, WWND, PRES, TPOT;

int i, j, k, iLevels, year, month, day, hour, min;

//    printf("Reading values...\n");
// printf("%s\n", infilename);

      if ( (infile = fopen(infilename, "r")) != NULL)
    {
        fgets(header, 300, infile);    // header data
       
        fgets(buffer, 120, infile);    // start date

        year = 2000 + atoi(&buffer[19]);
        month = atoi(&buffer[22]);
        day = atoi(&buffer[25]);
        hour = atoi(&buffer[28]);
        min = atoi(&buffer[31]);
        start_time = 10000*year+100*month+day+hour/100.+min/10000.;      // I admit this is very bad programming but I am returning the profile time and date as YYYYMMDD.HHMM

        fgets(buffer, 120, infile);    // end date
        fgets(buffer, 120, infile);    // blank
        fgets(buffer, 120, infile);    // line
        fgets(buffer, 120, infile);    // profile time

// Used Nearest Grid Point ( 245, 242) to Lat:    53.90, Lon:  -122.79
        fgets(buffer, 120, infile);    // nearest grid point
        lat = atof(&buffer[48]);
        lon = atof(&buffer[61]);

        fgets(buffer, 120, infile);    // 2D fields
        fgets(buffer, 120, infile);    // PRSS...
        fgets(buffer, 120, infile);    // hpa...

// Pressure is listed twice?
        fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &PRSS, &PRSS, &TPP1, &SHGT, &MXHT, &SHTF, &USTR, &DSWF);
// printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", PRSS, PRSS, TPP1, SHGT, MXHT, SHTF, USTR, DSWF);
        fgets(buffer, 120, infile);    // 3D fields

        fgets(buffer, 120, infile);    // UWND...

        fgets(buffer, 120, infile);    // m/s...


// UWND and VWND are listed twice.  Using the second set (W->E  S->N)
        if (strstr(header, "nest") !=NULL)
            fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                       &PRES, &UWND, &VWND, &TEMP, &SPHU, &WWND, &PRES, &TPOT, &UWND, &VWND);
        else
            fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                       &PRES, &PRES, &UWND, &VWND, & WWND, &TEMP, &SPHU, &TPOT, &UWND, &VWND); // first PRES ignored

// printf("... %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", PRES, UWND, VWND, TEMP, SPHU, WWND, PRES, TPOT, UWND, VWND);

        iLevels = 0;
        PREV_PRES = PRES;
        PREV_TEMP = TEMP;
         while( !feof( infile ) && iLevels < MAX_LEVELS )
        {
            if (iLevels == 0)
            {
                fbp->WS = 3.6 * sqrt(UWND*UWND + VWND*VWND);
                fbp->WD = 180./PI * atan2(UWND, VWND)+180.;
                feps->temp = TEMP;
                feps->rh=rh_calc(SPHU, PRES, TEMP);
                feps->ts = TEMP;
                feps->zs = SHGT; 
// printf("%lf %lf %lf %lf %lf %lf %lf %lf \n", UWND, VWND, fbp->WS, fbp->WD, feps->temp, feps->rh, feps->ts, feps->zs);
            }

            if (PRES < 850. && PREV_PRES >= 850.)
            {
                feps->t850 = PREV_TEMP + (TEMP-PREV_TEMP)/(PRES-PREV_PRES);             // temperature at 850 mb
                feps->z850 = 1457;            // height at 850 mb  - default value
            }

            if (PRES < 700. && PREV_PRES >= 700.)
            {
                feps->t700 = PREV_TEMP + (TEMP-PREV_TEMP)/(PRES-PREV_PRES);             // temperature at 700 mb
                feps->z700 = 3012;            // height at 700 mb  - default value
            }

            if (PRES < 500. && PREV_PRES >= 500.)
            {
                feps->t500 = PREV_TEMP + (TEMP-PREV_TEMP)/(PRES-PREV_PRES);             // temperature at 500 mb    
                feps->z500 = 5574;            // height at 500 mb  - default value
            }


            PREV_PRES = PRES;
            PREV_TEMP = TEMP;
            fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                            &PRES, &UWND, &VWND, &TEMP, &SPHU, &WWND, &PRES, &TPOT, &UWND, &VWND);

//            printf("%d %lf\n", iLevels, PRES);
            iLevels++;
        }
        
        fclose(infile);

        return (start_time);    // I admit this is very bad programming but I am returning the profile time and date as YYYYMMDD.HHMM
    }
    else
        return(-1);
}

//int Readupperair(char *infilename, FEPS *feps, double *P, double *Z, double *T, double *Td, double *W, double *D)
int Readupperair(char *infilename, FEPS *feps, UA *ua)
{
    int i,  ii, Hyb_Lvl, iLevels, type;
    double Pres, Hgt, TT, TD, Wspd, Wdir, phgt, ppres, ptt, ptd, UU, VV, RELH, MIXR, THTA, THTE, THTV,
          theta0, theta, zmix, sfcpres,
          zsfc, tsfc, tdsfc,
          z850, t850, td850,
          z700, t700, td700,
          z500, t500, td500;

    char buffer[128];

    
    if ((infile = fopen(infilename, "r")) == NULL)
    {
        printf("Error reading %s\n", infilename);
        iLevels = -999;
    }
    
    else
    {
/*
        ua->P = (double *)calloc(MAX_LEVELS, sizeof(double));
        ua->Z = (double *)calloc(MAX_LEVELS, sizeof(double));
        ua->T = (double *)calloc(MAX_LEVELS, sizeof(double));
        ua->Td = (double *)calloc(MAX_LEVELS, sizeof(double));
        ua->W = (double *)calloc(MAX_LEVELS, sizeof(double));
        ua->D = (double *)calloc(MAX_LEVELS, sizeof(double));     
*/

        ua->iLevels = 0;
        Hyb_Lvl = 0;                    
        ii = 0;
        type=0;   
        TT = 0;

        do
        {
            fgets(buffer, 640, infile);   
            ii++;
            
//            if (ii == 8 && strstr(buffer, "hPa")!=NULL)
            if (strstr(buffer, "hPa")!=NULL)
            {
                type = 8; // Univ of Wyoming: http://weather.uwyo.edu/upperair/sounding.html
                fgets(buffer, 640, infile); // -----------------------------------------------------------------------------
//                fscanf(infile, "%lf %lf", &Pres, &Hgt); //  1000.0     82                 
            }
            if (ii == 11 && strstr(buffer, "(mb)")!=NULL) type = 11; // Al's new profile format
            if (ii == 13 && strstr(buffer, "(mb)")!=NULL) type = 13; // Al's old profile format

        } while (type==0 && ii<13);

        if (type==8)
            do 
            {
                fgets(buffer, 640, infile);   
                if (strstr(buffer, "          ")==NULL)
                    sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Pres, &Hgt, &TT, &TD, &RELH, &MIXR, &Wdir, &Wspd, &THTA, &THTE, &THTV);
            } while (TT == 0);
        else if (type==11)
            fscanf(infile, "%d %lf %lf %lf %lf %lf %lf", &Hyb_Lvl, &Pres, &Hgt, &TT, &TD, &Wspd, &Wdir);   
        else
            fscanf(infile, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &Hyb_Lvl, &Pres, &Hgt, &TT, &TD, &UU, &VV, &Wspd, &Wdir);

        feps->zs = phgt = zsfc = Hgt;
        feps->pressure = sfcpres = ppres = Pres;
        feps->ts = ptt = tsfc = TT;
        ptd = tdsfc = TD;
        theta0 = (tsfc+273.16)*pow((1000./Pres), 0.286);
        zmix = 0;
  
        do
        {
            ua->P[ua->iLevels] = Pres;
            ua->Z[ua->iLevels] = Hgt;
            ua->T[ua->iLevels] = TT;
            ua->Td[ua->iLevels] = TD;
            ua->W[ua->iLevels] = Wspd;
            ua->D[ua->iLevels] = Wdir;

// printf("1... %d %lf %lf %lf %lf %lf %lf\n", Hyb_Lvl, Pres, Hgt, TT, TD, Wspd, Wdir);
// printf("2... %d %lf %lf %lf %lf %lf %lf\n", Hyb_Lvl, *P, *Z, *T, *Td, *W, *D);
            ua->iLevels++;

            if (Pres<850. && ppres>=850.)
            {
                feps->t850 = ptt + (TT-ptt)*(850-ppres)/(Pres-ppres);             // temperature at 850 mb
                feps->z850 = phgt + (Hgt-phgt)*(850-ppres)/(Pres-ppres);        // height at 850 mb
            }
            
            if (Pres<700. && ppres>=700.)
            {
                feps->t700 = ptt + (TT-ptt)*(700-ppres)/(Pres-ppres);             // temperature at 700 mb
                feps->z700 = phgt + (Hgt-phgt)*(700-ppres)/(Pres-ppres);        // height at 700 mb
            }
            
            if (Pres<500. && ppres>=500.)
            {
                feps->t500 = ptt + (TT-ptt)*(500-ppres)/(Pres-ppres);             // temperature at 500 mb
                feps->z500 = phgt + (Hgt-phgt)*(500-ppres)/(Pres-ppres);        // height at 500 mb
            }
        
            if (Pres<250. && ppres>=250.)
            {
                feps->t250 = ptt + (TT-ptt)*(500-ppres)/(Pres-ppres);             // temperature at 250 mb
                feps->z250 = phgt + (Hgt-phgt)*(500-ppres)/(Pres-ppres);        // height at 250 mb
            }
        
            theta = (TT+273.16)*pow((1000./Pres), 0.286);
        
            if (theta0 + 5 > theta) zmix = Hgt;
//          printf("%lf %lf %lf\n", theta0, theta, zmix);

            phgt = Hgt;
            ppres = Pres;
            ptt = TT;
            ptd = TD;


        if (type==8)
            fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Pres, &Hgt, &TT, &TD, &RELH, &MIXR, &Wdir, &Wspd, &THTA, &THTE, &THTV);
        else if (type==11)
            fscanf(infile, "%d %lf %lf %lf %lf %lf %lf", &Hyb_Lvl, &Pres, &Hgt, &TT, &TD, &Wspd, &Wdir);   
        else
            fscanf(infile, "%d %lf %lf %lf %lf %lf %lf %lf %lf", &Hyb_Lvl, &Pres, &Hgt, &TT, &TD, &UU, &VV, &Wspd, &Wdir);
        
        } while (Pres >= 100 && iLevels< MAX_LEVELS);
    
// printf("3... %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", sfcpres, feps->ts, feps->zs, feps->t850, feps->z850, feps->t700, feps->z700, feps->t500, feps->z500, zmix);    
// printf("4... %lf\n", (feps->t500 - feps->ts)/(feps->z500 - feps->zs));
        fclose(infile);
    }
    
    return(iLevels);
}

