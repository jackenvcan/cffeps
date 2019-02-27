#ifndef PFAS_H
#define PFAS_H

#define MAX_STNS 40
#define MAX_VALUES 100
#define ID_SIZE 24
#define NAME_SIZE 256
#define MAX_DMC 100
#define MAX_DAYS 100

typedef struct                        /* Weather station data*/
{
    char    *id[MAX_STNS],
            wxstnfilename[64];

    double  lat[MAX_STNS],
            lng[MAX_STNS],
            x[MAX_STNS],
            y[MAX_STNS],
            elev[MAX_STNS];

    int        wxstns;

} WXSTN;

typedef struct                        /* Weather station data*/
{
    char    *id[MAX_VALUES],
            valuefilename[64];

    double  lat[MAX_VALUES],
            lng[MAX_VALUES],
            x[MAX_VALUES],
            y[MAX_VALUES],
            elev[MAX_VALUES];

    int        values;

} VALUE;

typedef struct                        /* PFAS global data */
{
   double DMCex,
          idw_pow,
          friction,
          SlopeEffect,
          WxFactor,
          WnFactor,
          ws[MAX_STNS],
          pdf;

   int rad,
       timezone,
       model,
       smooth_days,
       smooth_dmcs;

   char diurnal[4];

}  PFAS;


typedef struct                        /* Risk data */
{
   int year,
       month,
       day,
       ndays,
       DMC,
       simdays[12];

   double minprob,
          friction; /* a universal friction value -- kinda obsolete now */

}  RISK;

typedef struct                        /* Simulation data */
{
   int year,
       month,
       day,
       ndays,
       DMC,
       direction,
       simdays[12],
       iterations;

   double lat,      /* latitude of ingition point */
          lng,      /* longitude of ingition point */
          minprob,  /* probability at wich the progagation routine ends */
          size,        /* use as an initial size of the fire */
          fsize,    /* final size: used for percentile calculations */
          trim,     /* m - used to trim the fuels data set */
          friction; /* a universal friction value -- kinda obsolete now */

   char igngrid[64],    /* path to an ASCII grid showing the ignition zone */
        hotspots[64],   /* path to a file of hotspots  */
        values[64],     /* a path to a file of values at risk  */
        method[12];     /* fire growth methodology: 16 point if method > 8, else 8 point,*/

}  SIM;

typedef struct
{
   int year,
    month,
    day,
    time,
    BE,
    accel,
    ensemble,
    slopeeffect,
    spotfire,
    firebreaks,
    extinction,
    direction,
    method,
        recursions,  
    simdays[12];

   char igngrid[64], 
        hotspots[64], /* a path to a file of hotspots */
        hfigrid[64], 
        cfbgrid[64], 
        sfcgrid[64], 
        tfcgrid[64];

   double lat,
          lng,
          ndays,
          FFMC,
          DMC,
          DC,
          BUI,
          GFL,
          C,
          CBH,
          PC,
          PDF,
          DMCex,
          timestep,
          size,        /* use as an initial size of the fire */
          trim,
          bias;

}  FGM;

typedef struct
{
   char diurnal[12],
        stream[12],
        file[64],
        forecast[64],
        beckcalc[64];

   int year,
       month,
       day,
       time,
       timezone,
       extrapolate,
       ndays,
       simdays[12],
       direction,
       ensemble,
       method,
       night,
       nrows, ncols;                /* if data is to be converted to a grid, nrows x ncols */

   double lat,
          lng,
          T,
          rh,
          ws,
          wd,
          Tn,
          Tx,
          RHn,
          Wn,
          Wx,
          WnFactor,
          WxFactor,
          FFMC,
          DMC,
          DC;

}  WX;

typedef struct
{
   double temp,
          rh,
          ws,
          wd,
          precip,
          FFMC,
          DMC,
          DC,
          ISI,
          BUI,
          FWI;
}  FWX;

typedef struct
{
   double temp[365],
          rh[365],
          ws[365],
          wd[365],
          precip[365],
          FFMC[365],
          DMC[365],
          DC[365],
          ISI[365],
          BUI[365],
          FWI[365];
}  FWXJ;

typedef struct
{
    char    *id[MAX_STNS];
    double    temp[MAX_STNS],
            rh[MAX_STNS],
            ws[MAX_STNS],
            wd[MAX_STNS],
            precip[MAX_STNS],
            FFMC[MAX_STNS],
            DMC[MAX_STNS],
            DC[MAX_STNS],
            ISI[MAX_STNS],
            BUI[MAX_STNS],
            FWI[MAX_STNS];
}  FWDAY;

typedef struct
{
   double temp[12][31],
          rh[12][31],
          ws[12][31],
          wd[12][31],
          precip[12][31],
          FFMC[12][31],
          DMC[12][31],
          DC[12][31],
          ISI[12][31],
          BUI[12][31],
          FWI[12][31];
}  FWYR;

typedef struct                        /* FWAP record */
{
   unsigned short dry,
                  wet,
                  humid,
                  wind,
                  rain,
                  dir,
                  ffmc,
                  dmc,
                  dc,
                  isi,
                  bui,
                  fwi,
                  dsr;
}  FWAP;

typedef struct
{
    double width[10], 
           HFImin[10];
}   FB;


typedef struct
{
   int DMCeqn[17];

   double depth[17],
          rhoB[17],
          Ri[17];

}  SMO;


typedef struct
{
    char method[12];

    int reps;

    double    t,
            rh,
            ws,
            wd;

}  ENS;

typedef struct
{
   char method[12], *name[100];

   int resources;

   double time,
          rate[100],
          width[100],
          cost[100],
          arrival[100],
          turnaround[100],
          area[100];

}  RES;



int Date2Julian(int, int, int);
void Julian2Date(int, int *, int *, int *);

void Error(char *);
double FindHFImin(double, FB *);

int ReadGlobalData(PFAS *ptr);
int ReadSimData(SIM *ptr);
int ReadFirebreaksData(FB *ptr);
int ReadFgmData(FGM *ptr);
int ReadResData(RES *ptr);
int ReadRiskData(RISK *ptr);
int ReadWxData(WX *ptr);
int ReadEnsembleData(ENS *ptr);
int ReadSmolderingData(SMO *ptr);
int PrintSimData(SIM *ptr);
int PrintFgmData(FGM *ptr);
int FBPType(char *, int, int, int);

/* System allows for 50 fuel types including non-fuels */
char *FBPs[50];
int iFuels[50], iCodes[50];
double *Psp[8], *Psp0[8];

static char *months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

static char *dirs[] = {"N", "NE", "E", "SE", "S", "SW", "W", "NW"};

static int jdays[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};                             

/* ASCII Grid Info */
int ncols, nrows, NODATA_value;
int C1,C2,C3,C4,C5,C6,C7,D1,M1,M2,M3,M4,O1a,O1b,S1,S2,S3,WA,NF, maxFuels;
double xllcorner, yllcorner, cellsize, missing;
int wxdays;
int iwidths;

/* filepath to data */
char filepath[64];
char errortext[64];
static double zero = 0.;
//static int TRUE=0, FALSE=-1;
#endif
