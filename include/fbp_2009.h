#ifndef FBP_2009_H
#define FBP_2009_H

/* FBP structure 
   information is passed to the FBP subroutine via this structure */

typedef struct
{
/* inputs */
   char   FuelType[4];
   int    Accel,           /* 1 = point, 0 = line (no longer accepts 2 as a line source) */
          Dj,              /* Julian Day */
          Do,              /* Julian day of minimum FMC */
          ELV,             /* Elevation [m ASL] */
          BUIEff;          /* BUI effect: 0 = yes, else no */
   double t,               /* Hours since ignition */
          FFMC,            /* FFMC */
          ISI,             /* ISI */
          BUI,             /* BUI */
          WS,              /* wind speed [kmh] */
          WD,              /* wind direction [degrees] */
          GS,              /* Slope [percent] */
          Aspect,          /* Aspect [degrees] */
          PC,              /* Percent Confier for M1/M2 */
          PDF,             /* Percent Dead Fir for M3/M4 */
          C,               /* Percent Cured for O1a/O1b (85% default) */
          GFL,             /* Grass Fuel Load [kg/m^2] (0.3 default) */
          CBH,             /* Crown to Base Height [m] (FBP defaults)*/
          CFL,             /* Crown Fuel Load [kg/m^2] (FBP defaults) */
          LAT,             /* Latitude [decimal degrees] */
          LON,             /* Longitude [decimal degrees] */
          FMC,             /* FMC if known */
          SH,              /* C6 Stand Height [m] - 2009 */
          SD,              /* C6 Stand Density [stems/ha] - 2009 */
          theta,           /* elliptical direction of calculation */

/* outputs */
          ROS,             /* Rate of Spread [m/min] */
          FROS,            /* Flank rate of Spread [m/min] */
          BROS,            /* Back Rate of Spread [m/min] */
          TROS,            /* Rate of Spread at angle theta [m/min] */
          HROSt,           /* Head Rate of Spread at time t [m/min] */
          FROSt,           /* Flank Rate of Spread at time t [m/min] */
          BROSt,           /* Back Rate of Spread at time t [m/min] */
          TROSt,           /* Rate of Spread at angle theta at time t [m/min] */
          CFB,             /* Crown Fraction Burned */
          FCFB,            /* Flank Crown Fraction Burned [%] */
          BCFB,            /* Back Crown Fraction Burned [%] */
          TCFB,            /* Crown Fraction Burned at angle thetea [%] */
          HFI,             /* Head Fire Intensity [kW/m] */
          FFI,             /* Head Fire Intensity [kW/m] */
          BFI,             /* Head Fire Intensity [kW/m] */
          TFI,             /* Head Fire Intensity [kW/m] */
          TFC,             /* Total Fuel Consumption [kg/m^2]  */
          FTFC,            /* Flank Total Fuel Consumption [kg/m^2]  */
          BTFC,            /* Back Total Fuel Consumption [kg/m^2]  */
          TTFC,            /* Total Fuel Consumption at angle theta [kg/m^2]  */
          SFC,             /* Surface Fuel Consumption [kg/m^2] */
          TI,              /* Time of Crown Fire initiation [hrs since ignition] */
          FTI,             /* Time of Flank Crown Fire initiation [hrs since ignition] */
          BTI,             /* Time of Back Crown Fire initiation [hrs since ignition] */
          TTI,             /* Time of Crown Fire initiation at angle theta [hrs since ignition] */
          LB,              /* Length to Breadth ratio */
          RAZ,             /* Spread direction azimuth */
          WSV;             /* Net vectored wind speed */
}  FBP;

/* Default values * /
	fbp.Accel=0;
	fbp.Dj=0;
	fbp.Do=0;
	fbp.ELV=0;
	fbp.BUIEff=0;
	fbp.t=0.;
	fbp.FFMC=0.;
	fbp.ISI=0.;
	fbp.BUI=0.;
	fbp.WS=0.;
	fbp.WD=0.;
	fbp.GS=0.;
	fbp.Aspect=0.;
	fbp.PC=0.;
	fbp.PDF=0.;
	fbp.C=85.;
	fbp.GFL=0.3;
	fbp.CBH=0.;
	fbp.LAT=0.;
	fbp.LON=0.;
	fbp.FMC=0.;
	fbp.theta=0.; */

static double PI = 3.1415926;
static double NoBUI = -1.;

/* Static values used by the FBP system 
   Note these values can be changed but must maintain internal consitency */
#define FUELTYPE_C1  0
#define FUELTYPE_C2  1
#define FUELTYPE_C3  2
#define FUELTYPE_C4  3
#define FUELTYPE_C5  4
#define FUELTYPE_C6  5
#define FUELTYPE_C7  6
#define FUELTYPE_D1  7
#define FUELTYPE_D2  8
#define FUELTYPE_M1  9
#define FUELTYPE_M2  10
#define FUELTYPE_M3  11
#define FUELTYPE_M4  12
#define FUELTYPE_S1  13
#define FUELTYPE_S2  14
#define FUELTYPE_S3  15
#define FUELTYPE_01A 16
#define FUELTYPE_01B 17
#define FUELTYPE_WA  18
#define FUELTYPE_NF  19

#define MAX_FUELS    20

/* Note that for the 2009 version a, b and c terms ave been added for M3 and M4, though the RSI calculation involves D1 fule type as well in the revised eqns 29 to 32 -  KRA 12/13/2010 */
/* for ECCC application:
   Jack add D2(parameters)=D1 - D2 is unofficial type being D1 after leafs out also set WA and NF parameters to 0 */
static char* Fuels[MAX_FUELS]=
{  "C1",  "C2",   "C3",  "C4",  "C5",  "C6",  "C7",  "D1",  "D2",  "M1",  "M2",  "M3",  "M4",  "S1",  "S2",  "S3", "O1a",  "O1b", "WA", "NF"};
static double a[MAX_FUELS] =
{   90.,  110.,   110.,  110.,   30.,   30.,   45.,   30.,   30.,     0,     0,  120.,  100.,   75.,   40.,   55.,  190.,   250.,    0,   0 };
static double b[MAX_FUELS] =
{ .0649, .0282,  .0444, .0293, .0697, .0800, .0305, .0232, .0232,     0,     0, .0572, .0404, .0297, .0438, .0829, .0310,  .0350,    0,   0 };
static double c[MAX_FUELS] =
{   4.5,   1.5,    3.0,   1.5,   4.0,   3.0,   2.0,   1.6,   1.6,     0,     0,   1.4,  1.48,   1.3,   1.7,   3.2,   1.4,    1.7,    0,   0 };
static double q[MAX_FUELS] =
{  0.90,  0.70,   0.75,  0.80,  0.80,  0.80,  0.85,  0.90,  0.90,    .8,    .8,    .8,    .8,  0.75,  0.75,  0.75,  1.00,   1.00,    0,   0 };
static double BUIo[MAX_FUELS] =
{   72.,   64.,    62.,   66.,   56.,   62.,  106.,   32.,   32.,    50,    50,    50,    50,   38.,   63.,   31.,    1.,     1.,    0,   0 };
static double CBHs[MAX_FUELS] =
{    2.,    3.,     8.,    4.,   18.,    7.,   10.,     0,     0,    6.,    6.,    6.,    6.,     0,     0,     0,     0,      0,    0,   0 };
static double CFLs[MAX_FUELS] =
{   .75,   .80,   1.15,  1.20,  1.20,  1.80,   .50,     0,     0,    .8,    .8,    .8,    .8,     0,     0,     0,     0,      0,    0,   0 };
/* These values taken from Anderson 2002 */
static double depths[MAX_FUELS] =
{   3.4,  10.0,    6.5,   6.2,   4.6,   5.0,   5.0,   2.4,   2.4,   5.0,   5.0,   7.5,   7.5,   7.4,   7.4,   7.4,     0,      0,    0,   0 };
static double rhoBs[MAX_FUELS] =
{ 0.045, 0.034,   0.02, 0.031, 0.093,  0.05,  0.02, 0.061, 0.061, 0.108, 0.108, 0.061, 0.061, 0.078, 0.132,   0.1,     0,      0,    0,   0 };
static double Ris[MAX_FUELS] =
{  0.05,   0.0,   0.15,  0.15,  0.15,  0.15,  0.15,  0.59,  0.59,  0.25,  0.25,  0.15,  0.15,  0.15,  0.15,  0.15,     0,      0,    0,   0 };





/************************ fbp97.c *****************************/
/*    FBPCalc conducts all the calculations required for the 
      Canadian Forest Fire Behavior Prediction (FBP) system.
      The FBP system is documented in 

      Forestry Canada Fire Danger Group.  1992.  Development and
         Structure of the Canadian Forest Fire Behavior Prediction
         System.  For. Can., Sci. Sustainable Develop. Directorate, 
         Ottawa, Ont, Inf. Rep. ST-X-3.  63 p.

      Variable names and equation numbers are consistent with those 
      used in the FBP document.

      Note that FBPCalc is the principle call to FBP.  Inputs are 
      loaded into the FBP structure and FBPCalc is called. Outputs 
      are stored in the FBP stucture.

      The other subroutines are called within FBPCalc but are 
      listed here for completeness */

int FBPCalc(FBP *ptr);

/* FBPFuel determines an index for the FBP fuel type used throughout the program 
   Note that the order can be adjusted by varying the order in the static variables */
   int FBPFuel(char * FuelType);

/* Foliar Moisture Content (FMC) calculation 
   Note that 0.5 is added before the integer conversion in equations 2 and 4
   Note that equations 1 and 3 use positive longitude values for Canada */
   double FMCcalc(double LAT, double LON, int ELV, int Dj, int Do);

/* Surface Fuel Consumption (SFC) calculation */
   double SFCcalc(int Fuel, double FFMC, double BUI, double PC, double GFL);

/* Rate of Spread calculations */
   double ROScalc(int Fuel, double ISI, double BUI, double FMC, double SFC, 
                  double PC, double PDF, double C, double CBH);

/* Initial Spread Index (ISI) calculations -- used in slope effect calucations */
   double ISICalc(double FFMC, double WSV);


/* Effect of Slope on Rate of Spread */
   void Slopecalc(int Fuel, double FFMC, double BUI, double WS, double WAZ, double GS, 
                  double SAZ, double FMC, double SFC, double PC, double PDF, double C, double CBH,
                  double *pRAZ, double *pWSV);

/* BUI Effect */
   double BEcalc(int Fuel, double BUI);

/* Crown Fraction Burned (CFB) calculation */
   double CFBcalc(int Fuel, double FMC, double SFC, double ROS, double CBH);

/* C6 has its own crowning fire model; hence, it has a special subroutine */
   void C6calc(int Fuel, double ISI, double BUI, double FMC, double SFC, double CBH, double *ROS, double *CFB, double *RSC);

/* Total Fuel Consumption (TFC) calculation */
   double TFCcalc(int Fuel, double CFL, double CFB, double SFC, double PC, double PDF);

/* Fire Intensity (*FI) calculation */
   double FIcalc(double FC, double ROS);

/* Back Rate of Spread (BROS) calculation */
   double BROScalc (int Fuel, double FFMC, double BUI, double WSV, double FMC, 
                    double SFC, double PC, double PDF, double C, double CBH);

/* Rate of Spread at time t (since ignition) calculation */
   double ROStcalc(int Fuel, double ROSeq, double t, double CFB);

/* Length to Breadth ratio (LB) calculation */
   double LBcalc(int Fuel, double WSV);

/* Length to Breadth ratio (LB) calculation  with time and acceleration */
   double LBtcalc(int Fuel, double LB, double t, double CFB);

/* Flank Rate of Spread (FROS) calculation */
   double FROScalc(double ROS, double BROS, double LB);

/* Rate of Spread (TROS) calculation with respect to departure from wind direction */
   double ROSthetacalc(double ROS, double FROS, double BROS, double theta);

/* FBPdefaults loads default values into the FBP structure (mostly zeros).
   This is useful to avoid unexpected errors due to oversights */
   int FBPdefaults(FBP *fbp);

/* FBP test runs this code agains the expected results as described by BMW */
   void FBPtest();

/************************* fbp2009pro.c ***************************************/

/*    From one non-zero fire shape input parameter (area, perimeter, forward 
      spread distance, time1, time2), PROCalc calculates all other fire shape
      parameters.  Accelleration effects may be included.

      Note that all desired outputs must be filled with 0 when called. If all
      input parameters are 0, the outputs are zero.

      If the headfire rate of spread (Rh) is zero, then the second elapsed time
      is set to zero and all remaining parameters are untouched (presumably
      all but one have already been set to zero by the calling routine).  */

void PROCalc(char   *Fuel,       /* The only significant fuel type is 4 (C1)        */
             int    Accel,       /* Are Acceleration effects be incorporated? (T/F) */
             double CFB,         /* Crown Fraction Burned                           */
             double Rh,          /* Headfire Eq. Rate of spread in m/min            */
             double Rf,          /* Flankfire Eq. Rate of spread in m/min           */
             double Rb,          /* Backfire Eq. Rate of spread in m/min            */
             double *t1,         /* Elapsed time 1 in hrs                           */
             double *t2,         /* Elapsed time 2 in hrs                           */
             double *A,          /* Area burned in heactares                        */
             double *P,          /* Perimeter encompassed in kilometres             */
             double *D );        /* Distance travelled in kilometres                */

/************************* THE END ***************************************/

#endif
