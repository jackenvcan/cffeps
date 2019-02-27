/*
 * fbp_2009.c
 *
 *  Created on: Jan 6, 2011
 *      Author: Kerry Anderson
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fbp_2009.h"

// int Date2Julian();

FILE *infile, *outfile5, *outfile6, *intfile;

/********************** Kerry Anderson's Stuff ***************************/

/* This subroutine represents equations from

Wotton, B.M.; Alexander, M.E.; Taylor, S.W.  2009,.  Updates and revisions to the 
    1992 Canadian Forest Fire Behavior Prediction System .  Natural Resources 
    Canada, Canadian Forest Service, Great Lakes Forestry Centre, Sault Ste. Marie, 
    Ontario, Canada.  Infomration Report GLC-X-10, 45 p.

Updates are indicated with a "- 2009" in the equation number comment on the right

*/

/* FBP Structure is stored in fbp09.h

   Typical call from C:  

   FBP fbp;   (declaration of the fbp variable structure)
   status = FBPCalc(&fbp);

   Discrepancies between BMW's fbp.c and ST-X-3 

   eqn 32: BMW has 33.5 instead of 35.5 (see ROScalc and Slopecalc)
   eqn 57: BMW had forced SFC to 2.0 for C6 (this is now gone)
   b[O1b]: BMW has 0.0829 instead of 0.0310
   c[O1a]: BMW has 1.41 instead of 1.4
*/

// Note: There is no attempt to account for greened-up deciduous (D2)

int FBPCalc (FBP *fbp)
{

/*   printf("%s\n%d %d %d %d %d\n%lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf\n%lf %lf %lf %lf %lf\n%lf %lf\n", 
        fbp->FuelType, 

        fbp->Accel,
        fbp->Dj,
        fbp->Do,
        fbp->ELV,
        fbp->BUIEff,

        fbp->t,
        fbp->FFMC,
        fbp->ISI,
        fbp->BUI,
        fbp->WS,
        
        fbp->WD,
        fbp->GS,
        fbp->Aspect,
        fbp->PC,
        fbp->PDF,
        
        fbp->C,
        fbp->GFL,
        fbp->CBH,
        fbp->LAT,
        fbp->LON,
        
        fbp->FMC,
        fbp->theta); */
 
/* inputs */
   char   FuelType[4]; 
   int    Accel,           /* 0 = point, 1 = line (no longer accepts 2 as a line source) */
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

/* intermediate variables */
double E, LBt, ROSt, RSC, SAZ, WAZ;

int Fuel, status;

    status = 0;    // success
        
   strcpy(FuelType, fbp->FuelType); 

   Accel  = fbp->Accel;
   Dj     = fbp->Dj;
   Do     = fbp->Do;
   ELV    = fbp->ELV;
   BUIEff = fbp->BUIEff;
   t      = fbp->t;
   FFMC   = fbp->FFMC;
   ISI    = fbp->ISI;
   BUI    = fbp->BUI;
   WS     = fbp->WS;
   WD     = fbp->WD * PI/180.;        /* radians */
   GS     = fbp->GS;
   Aspect = fbp->Aspect * PI/180.;    /* radians */
   PC     = fbp->PC;
   PDF    = fbp->PDF;
   C      = fbp->C;
   GFL    = fbp->GFL;
   CBH    = fbp->CBH;
   CFL    = fbp->CFL;
   LAT    = fbp->LAT;
   LON    = fbp->LON;
   FMC    = fbp->FMC;
   theta  = fbp->theta * PI/180.;    /* radians */
   SD     = fbp->SD;
   SH     = fbp->SH;
   
/*   
printf("%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",  
        Accel, Dj, Do, ELV, BUIEff, 
        t, FFMC, ISI, BUI, WS, WD, GS, Aspect, PC, PDF, C, GFL, CBH, CFL, LAT, LON, FMC, theta);
*/
/* Some default values */
      if (Accel <= 0) Accel=0; else Accel=1;    /* line=1 (no accelleration effect) */    /* changed from < 0 to <= 0 2017-03-20 */
      if (Dj < 0 || Dj > 366) Dj = 0;
      if (Do < 0 || Do > 366) Do = 0;
      if (ELV < 0.0 || ELV > 10000.) ELV = 0.;
      if (BUIEff < 0.0) BUIEff=-1; else BUIEff=0;
      if (t < 0.0) t = -t;
      if (t > 366*24.) t = 24.;
      if (FFMC < 0.0 || FFMC > 101.) FFMC = 0.0;
      if (ISI < 0.0 || ISI > 300.) ISI = 0.0;
      if (BUI < 0.0 || BUI > 1000.) BUI = 0.0;
      if (WS < 0.0 || WS > 300.) WS = 0.0;
      if (WD < -2*PI || WD > 2*PI) WD = 0.0;
      if (GS < 0.0 || GS > 200.) GS = 0.0;
      if (Aspect < -2*PI || Aspect > 2*PI) GS = 0.0;
      if (PC < 0.0 || PC > 100.) PC = 50.0;
      if (PDF < 0.0 || PDF > 100.) PDF = 35.0;
      if (C <= 0.0 || C > 100.) C = 95.;
      if (GFL <= 0.0 || GFL > 100.) GFL = 0.35;  /* changed from 0.3 to 0.35, pg 6 - 2009 */
      if (LAT < -90.0 || LAT > 90.) LAT = 0.0;
      if (LON < -180. || LON > 360.) LON = 0.0;
      if (theta < -2*PI || theta > 2*PI) theta = 0.0;
      if (SD <0. || SD > 100000.) SD=-999.;
      if (SH <0. || SH > 100.) SH=-999.;
      

/* Convert time from hours to minutes */
   t = t*60.;

/* Corrections to reorient WAZ, SAZ */
   WAZ = WD + PI;
   if (WAZ > 2*PI) WAZ = WAZ - 2*PI;

/* nb: BMW's data set appears to have aspect not saz */  
   SAZ = Aspect + PI;
   if (SAZ > 2*PI) SAZ = SAZ - 2*PI;

/* Make LON positive for the Western Hemisphere */
   LON = -LON;

/* Enter FBP Calculations */
   Fuel = FBPFuel(FuelType);
   if(Fuel==-1 || FFMC == 0. || BUI == 0.)
   {
      fbp->ROS   = 0.;
      fbp->FROS  = 0.;
      fbp->BROS  = 0.;
      fbp->TROS  = 0.;
      fbp->HROSt = 0.;
      fbp->FROSt = 0.;
      fbp->BROSt = 0.;
      fbp->TROSt = 0.;
      fbp->CFB   = 0.;
      fbp->FCFB  = 0.;
      fbp->BCFB  = 0.;
      fbp->TCFB  = 0.;
      fbp->HFI   = 0.;
      fbp->FFI   = 0.;
      fbp->BFI   = 0.;
      fbp->TFI   = 0.;
      fbp->TFC   = 0.;
      fbp->FTFC  = 0.;
      fbp->BTFC  = 0.;
      fbp->TTFC  = 0.;
      fbp->SFC   = 0.;
      fbp->TI    = -999.;
      fbp->FTI   = -999.;
      fbp->BTI   = -999.;
      fbp->TTI   = -999.;
      fbp->LB    = -999.;
      fbp->RAZ    = -999.;
      fbp->WSV    = -999.;
      status = -1;    // fail

      return(status);
   }

   CFB=0.;


/* presently, we do not accept a zero CBH; use near zero if necessary */
   if (CBH <= 0. || CBH > 50.) 
      if (!strcmp(FuelType, "C6") && SD > 0. && SH > 0.) 
      {
         CBH = -11.2 + 1.06*SH + 0.00170*SD;                        /* 91 */
         if (CBH <0.) CBH = 0.;
      }
      else
         CBH = CBHs[Fuel];

/* presently, we do not accept a zero CFL,; use near zero if necessary */
   if (CFL <= 0.0 || CFL > 2.0) 
      CFL = CFLs[Fuel];

   if (FMC <= 0 || FMC > 120.)
      FMC = FMCcalc(LAT, LON, ELV, Dj, Do);

   SFC = SFCcalc(Fuel, FFMC, BUI, PC, GFL);

   if (!BUIEff) BUI = 0.0;  /* This turns off BUI effect */

   if (GS > 0. && FFMC > 0.)
   {
      Slopecalc(Fuel, FFMC, BUI, WS, WAZ, GS, 
                SAZ, FMC, SFC, PC, PDF, C, CBH,  
                &RAZ, &WSV);
   }
   else
   {
      WSV = WS;
      RAZ = WAZ;
   }

   if (FFMC > 0.)
      ISI = ISICalc(FFMC, WSV);

 if (!strcasecmp(Fuels[Fuel], "C6"))
      C6calc(Fuel, ISI, BUI, FMC, SFC, CBH,  &ROS, &CFB, &RSC);/* We use C6calc to calculate CFB */
   else
   {

      ROS = ROScalc(Fuel, ISI, BUI, FMC, SFC, PC, PDF, C, CBH);
      if (CFL > 0.)
         CFB  = CFBcalc(Fuel, FMC, SFC, ROS, CBH);
      else
         CFB = 0.;
   }


 LB = LBcalc(Fuel, WSV);

   if(Accel)
      LBt = LB;
   else
      LBt = LBtcalc(Fuel, LB, t, CFB);

   BROS = BROScalc(Fuel, FFMC, BUI, WSV, FMC, SFC, PC, PDF, C, CBH);
   FROS = FROScalc(ROS, BROS, LB);

/* TROS is the rate of spread towards angle theta */
   E = sqrt(1.-1./LB/LB);      /* eccentricity */
   TROS= ROS * (1.-E)/
          (1.-E*cos(theta - RAZ)); /* note: this is the old method using the focus as the ignition point */

//   TROS = ROSthetacalc(ROS, FROS, BROS, theta);

   if(Accel)
   {
      ROSt  = ROS;
      FROSt = FROS;
      BROSt = BROS;
      TROSt = TROS;
   }
   else
   {
      ROSt  = ROStcalc(Fuel, ROS, t, CFB);
      BROSt = ROStcalc(Fuel, BROS, t, CFB);
      FROSt = FROScalc(ROSt, BROSt, LBt);
//      TROSt = ROSthetacalc(ROSt, FROSt, BROSt, theta);
      
   E = sqrt(1.-1./LBt/LBt);      /* eccentricity */
   TROSt= ROSt * (1.-E)/
          (1.-E*cos(theta - RAZ)); /* note: this is the old method using the focus as the ignition point */
   }

   if (CFL == 0.)
   {
      CFB = 0.;
      FCFB = 0.;
      BCFB = 0.;
      TCFB = 0.;
   }
   else
   {
      if (!strcasecmp(Fuels[Fuel], "C6"))  
      {
         FCFB = 0.;
         BCFB = 0.;
         TCFB = 0.;
      }
      else
      {
/* equilibrium values */
         CFB  = CFBcalc(Fuel, FMC, SFC, ROS, CBH);
         FCFB = CFBcalc(Fuel, FMC, SFC, FROS, CBH);
         BCFB = CFBcalc(Fuel, FMC, SFC, BROS, CBH);
         TCFB = CFBcalc(Fuel, FMC, SFC, TROS, CBH);
      }
   }

   TFC  = TFCcalc(Fuel, CFL, CFB, SFC, PC, PDF);
   FTFC = TFCcalc(Fuel, CFL, FCFB, SFC, PC, PDF);
   BTFC = TFCcalc(Fuel, CFL, BCFB, SFC, PC, PDF);
   TTFC = TFCcalc(Fuel, CFL, TCFB, SFC, PC, PDF);

/* equilibrium values */
   HFI = FIcalc(TFC,  ROS);
   FFI = FIcalc(FTFC, FROS);
   BFI = FIcalc(BTFC, BROS);
   TFI = FIcalc(TTFC, TROS);

/* For now... */
   TI = 0.;
   FTI = 0.;
   BTI = 0.;
   TTI = 0.;

   fbp->ROS   = ROS;
   fbp->FROS  = FROS;
   fbp->BROS  = BROS;
   fbp->TROS  = TROS;

   if (fbp->t < 0)
   {
        fbp->HROSt = -ROSt;
        fbp->FROSt = -FROSt;
        fbp->BROSt = -BROSt;
        fbp->TROSt = -TROSt;
        fbp->CFB   = -CFB;
   }
   else
   {
        fbp->HROSt = ROSt;
        fbp->FROSt = FROSt;
        fbp->BROSt = BROSt;
        fbp->TROSt = TROSt;
   }

   fbp->CFB  = CFB;
   fbp->FCFB  = FCFB;
   fbp->BCFB  = BCFB;
   fbp->TCFB  = TCFB;
   fbp->HFI   = HFI;
   fbp->FFI   = FFI;
   fbp->BFI   = BFI;
   fbp->TFI   = TFI;
   fbp->TFC   = TFC;
   fbp->FTFC  = FTFC;
   fbp->BTFC  = BTFC;
   fbp->TTFC  = TTFC;
   fbp->SFC   = SFC;
   fbp->TI    = TI;
   fbp->FTI   = FTI;
   fbp->BTI   = BTI;
   fbp->TTI   = TTI;
   fbp->LB    = LB;
   fbp->RAZ   = RAZ * 180./PI;
   fbp->WSV   = WSV;
   fbp->FMC   = FMC;    /* failed to return FMC: 5/5/2014 */

   return(status);
}

/* Determine and index for the FBP fuel type used throughout the program 
   Note that the order can be adjusted by varying the order in the static variables */

int FBPFuel(char * FuelType)
{
int i, Fuel;
char temp[4];

   if (!strcasecmp(FuelType, "O1"))
      strcpy(temp, "O1b");
   else
      strcpy(temp, FuelType);

   Fuel = -1;    // fail
   for (i = 0; i < MAX_FUELS; i++)
      if (!strcasecmp(FuelType, Fuels[i]) && strcasecmp(FuelType, "WA") && strcasecmp(FuelType, "NF")) Fuel = i;

   return(Fuel);
}

/* Foliar Moisture Content (FMC) calculation 
   Note that 0.5 is added before the integer conversion in equations 2 and 4
   Note that equations 1 and 3 use positive longitude values for Canada */

double FMCcalc(double LAT, double LON, int ELV, int Dj, int Do)
  {
    double FMC, LATN;
    int ND;

    FMC = -1.;

   if (Do <= 0)      /* if Do, date of min FMC, is not known then Do <= 0 */
   {
      if (ELV <= 0)
      {
         LATN = 46.0 + 23.4 * exp(-0.0360 *(150. - LON));              /* 1 */
         Do = (int)(151. * LAT/LATN + 0.5);                      /* 2 (+0.5) */
      }
      else
      {
         LATN = 43. + 33.7 * exp(-0.0351*(150 - LON));                 /* 3 */
         Do = (int)(142.1 * LAT/LATN + (0.0172 * ELV) + 0.5);    /* 4 (+0.5) */
      }
   }

   ND = abs(Dj - Do);                                                  /* 5 */
   
   if (ND < 30)
      FMC = 85. + 0.0189 * ND * ND;                                    /* 6 */
   else 
      if (ND < 50)
         FMC = 32.9 + 3.17 * ND - 0.0288 * ND * ND;                    /* 7 */
      else
         FMC  = 120.;                                                  /* 8 */

   return(FMC);
}

/* Surface Fuel Consumption (SFC) calculation */

double SFCcalc(int Fuel, double FFMC, double BUI, double PC, double GFL)
{
double SFC, FFC, WFC;

   SFC = -1.;
   FFC = 0.;
   WFC = 0.;

   if (!strcasecmp(Fuels[Fuel], "C1"))
//      SFC = 1.5 * (1. - exp(-0.230 * (FFMC - 81)));                    /* 9 */
      if (FFMC > 84.)
         SFC = 0.75 + 0.75*pow((1.-exp(-0.23*(FFMC-84.))), 0.5);       /* 9a - 2009 */
      else
         SFC = 0.75 - 0.75*pow((1.-exp(-0.23*(84.-FFMC))), 0.5);       /* 9b - 2009 */


   if (!strcasecmp(Fuels[Fuel], "C2") ||
       !strcasecmp(Fuels[Fuel], "M3") ||
       !strcasecmp(Fuels[Fuel], "M4"))
      SFC = 5.0 * (1. - exp(-0.0115 * BUI));                           /* 10 */

   if (!strcasecmp(Fuels[Fuel], "C3") ||
       !strcasecmp(Fuels[Fuel], "C4"))
      SFC = 5.0 * pow(1. - exp(-0.0164*BUI), 2.24);                    /* 11 */

   if (!strcasecmp(Fuels[Fuel], "C5") ||
       !strcasecmp(Fuels[Fuel], "C6"))
      SFC = 5.0 * pow(1. - exp(-0.0149*BUI), 2.48);                    /* 12 */

   if (!strcasecmp(Fuels[Fuel], "C7"))
   {
      if (FFMC> 70.)
         FFC = 2. * (1. - exp(-0.104 * (FFMC - 70)));                  /* 13 */
      else
         FFC = 0.;
      WFC = 1.5 * (1. - exp(-0.0201 * BUI));                           /* 14 */
      SFC = FFC + WFC;                                                 /* 15 */
   }

   if (!strcasecmp(Fuels[Fuel], "D1"))
      SFC = 1.5 * (1. - exp(-0.0183 * BUI));                           /* 16 */

   if (!strcasecmp(Fuels[Fuel], "M1") ||
       !strcasecmp(Fuels[Fuel], "M2"))
   {
      SFC = (PC/100. * SFCcalc(FBPFuel("C2"), FFMC, BUI, PC, GFL))
      + ((100. - PC)/100*SFCcalc(FBPFuel("D1"), FFMC, BUI, PC, GFL));  /* 17 */
   }

   if (!strcasecmp(Fuels[Fuel], "O1a") ||
       !strcasecmp(Fuels[Fuel], "O1b"))
      SFC = GFL;                                                       /* 18 */

   if (!strcasecmp(Fuels[Fuel], "S1"))
   {
      FFC = 4.0 * (1. - exp(-0.025*BUI));                              /* 19 */
      WFC = 4.0 * (1. - exp(-0.034*BUI));                              /* 20 */
      SFC = FFC + WFC;                                                 /* 25 */
   }

   if (!strcasecmp(Fuels[Fuel], "S2"))
   {
      FFC = 10.0 * (1. - exp(-0.013*BUI));                             /* 19 */
      WFC = 6.0 * (1. - exp(-0.060*BUI));                              /* 20 */
      SFC = FFC + WFC;                                                 /* 25 */
   }

   if (!strcasecmp(Fuels[Fuel], "S3"))
   {
      FFC = 12.0 * (1. - exp(-0.0166*BUI));                            /* 19 */
      WFC = 20.0 * (1. - exp(-0.0210*BUI));                            /* 20 */
      SFC = FFC + WFC;                                                 /* 25 */
   }

   if (SFC <= 0.) SFC = 0.000001;

   return(SFC);
}

/* Rate of Spread calculations */

double ROScalc(int Fuel, double ISI, double BUI, double FMC, double SFC, double PC, double PDF, double C, double CBH)
{
double RSI, RSI_M3, RSI_M4, CF, BE, ROS, CFB, RSC;


   RSI = -1.;

/* Note that only preliminary RSS calculations are done for C6 in this routine */
   if (!strcasecmp(Fuels[Fuel], "C1") || 
       !strcasecmp(Fuels[Fuel], "C2") ||
       !strcasecmp(Fuels[Fuel], "C3") ||
       !strcasecmp(Fuels[Fuel], "C4") ||
       !strcasecmp(Fuels[Fuel], "C5") ||
       !strcasecmp(Fuels[Fuel], "C7") ||
       !strcasecmp(Fuels[Fuel], "D1") ||
       !strcasecmp(Fuels[Fuel], "S1") ||  
       !strcasecmp(Fuels[Fuel], "S2") ||
       !strcasecmp(Fuels[Fuel], "S3")  )
      RSI = a[Fuel] * pow(1. - exp(-b[Fuel] * ISI), c[Fuel]);          /* 26 */

   if (!strcasecmp(Fuels[Fuel], "M1"))
      RSI =       PC/100. * ROScalc(FBPFuel("C2"), ISI, NoBUI, FMC, SFC, PC, PDF, C, CBH)
         + (100.-PC)/100. * ROScalc(FBPFuel("D1"), ISI, NoBUI, FMC, SFC, PC, PDF, C, CBH);  /* 27 */

   if (!strcasecmp(Fuels[Fuel], "M2"))
      RSI =       PC/100. * ROScalc(FBPFuel("C2"), ISI, NoBUI, FMC, SFC, PC, PDF, C, CBH)
      +0.2*(100.-PC)/100. * ROScalc(FBPFuel("D1"), ISI, NoBUI, FMC, SFC, PC, PDF, C, CBH);  /* 28 */

   if (!strcasecmp(Fuels[Fuel], "M3"))
   {
      RSI_M3 = a[Fuel] * pow(1. - exp(-b[Fuel] * ISI), c[Fuel]);       /* 30 - 2009 */

      RSI = PDF/100.* RSI_M3 + 
           (1.-PDF/100.)* ROScalc(FBPFuel("D1"), ISI, NoBUI, FMC, SFC, PC, PDF, C, CBH);    /* 29 - 2009 */
   }

   if (!strcasecmp(Fuels[Fuel], "M4"))
   {
      RSI_M4 = a[Fuel] * pow(1. - exp(-b[Fuel] * ISI), c[Fuel]);       /* 32 - 2009 */

      RSI = PDF/100.* RSI_M4 + 
      0.2*(1.-PDF/100.)* ROScalc(FBPFuel("D1"), ISI, NoBUI, FMC, SFC, PC, PDF, C, CBH);     /* 31 - 2009 */
   }

   if (!strcasecmp(Fuels[Fuel], "O1a") ||
       !strcasecmp(Fuels[Fuel], "O1b"))
   {

         if (C < 58.8)
            CF = 0.005*(exp(0.061*C)-1.);                              /* 35a - 2009 */
         else
            CF = 0.176 + 0.02*(C-58.8);                                /* 35b - 2009 */

/* RSI has been substituted for ROS in eqn 36 */
         RSI = a[Fuel] * pow(1. - exp(-b[Fuel] * ISI), c[Fuel]) * CF;  /* 36 */
   }

   if (!strcasecmp(Fuels[Fuel], "C6"))
      C6calc(Fuel, ISI, BUI, FMC, SFC, CBH,  &ROS, &CFB, &RSC);/* included here for completeness */
   else
   {
      BE = BEcalc(Fuel, BUI);
      ROS = BE*RSI;
   }

   if (ROS <= 0.) ROS = 0.000001;

   return(ROS);
}

double ISICalc(double FFMC, double WSV)

{
  double fW,m,fF,ISI;
     
     m = 147.2*(101.-FFMC)/(59.5+FFMC);                                /* 46 */
     fF = 91.9*exp(-0.1386*m)*(1.+pow(m,5.31)/4.93e7);                 /* 45 */

   if (WSV < 40.)
      fW = exp(0.05039*WSV);                                           /* 53 */
   else
      fW = 12. * (1. - exp(-0.0818 * (WSV-28.)));                      /* 53a */

     ISI = 0.208*fW*fF;                                                /* 52 */

  return(ISI);
}

/* Effect of Slope on Rate of Spread */

void Slopecalc(int Fuel, double FFMC, double BUI, double WS, double WAZ, double GS, 
               double SAZ, double FMC, double SFC, double PC, double PDF, double C, double CBH,
               double *pRAZ, double *pWSV)
{
double SF, RSZ, RSF, RSF_C2, RSF_D1, RSF_M3, RSF_M4, PDF100, ISZ, ISF, ISF_C2, ISF_D1, ISF_M3, ISF_M4, CF, m, fF, WSE, WSX, WSY; 

   if (GS >= 70.) 
      SF=10.;
   else
      SF = exp(3.533 * pow(GS/100., 1.2));                             /* 39 */

   ISZ = ISICalc(FFMC, 0.);
   RSZ = ROScalc(Fuel, ISZ, NoBUI, FMC, SFC, PC, PDF, C, CBH);
   RSF = RSZ * SF;                                                     /* 40 */

   if (!strcasecmp(Fuels[Fuel], "C1") || 
       !strcasecmp(Fuels[Fuel], "C2") ||
       !strcasecmp(Fuels[Fuel], "C3") ||
       !strcasecmp(Fuels[Fuel], "C4") ||
       !strcasecmp(Fuels[Fuel], "C5") ||
       !strcasecmp(Fuels[Fuel], "C6") ||
       !strcasecmp(Fuels[Fuel], "C7") ||
       !strcasecmp(Fuels[Fuel], "D1") ||
       !strcasecmp(Fuels[Fuel], "S1") ||
       !strcasecmp(Fuels[Fuel], "S2") ||
       !strcasecmp(Fuels[Fuel], "S3"))

   if ((1. - pow(RSF/a[Fuel], 1./c[Fuel])) >= 0.01)
      ISF = log(1. - pow(RSF/a[Fuel], 1./c[Fuel]))/(-b[Fuel]);      /* 41a - 2009 */
   else
      ISF = log(.01)/(-b[Fuel]);                                    /* 41b - 2009 */

   if (!strcasecmp(Fuels[Fuel], "M1") ||
       !strcasecmp(Fuels[Fuel], "M2"))
   {
      RSZ = ROScalc(FBPFuel("C2"), ISZ, NoBUI, FMC, SFC, PC, PDF, C, CBH);
      RSF_C2 = RSZ * SF;                                                     /* 40 */
      RSZ = ROScalc(FBPFuel("D1"), ISZ, NoBUI, FMC, SFC, PC, PDF, C, CBH);
      RSF_D1 = RSZ * SF;                                                     /* 40 */

      if ((1. - pow(RSF_C2/a[FBPFuel("C2")], 1./c[FBPFuel("C2")])) >= 0.01)
         ISF_C2 = log(1. - pow(RSF_C2/a[FBPFuel("C2")], 1./c[FBPFuel("C2")]))
               /(-b[FBPFuel("C2")]);                                   /* 41a - 2009 */
      else
         ISF_C2 = log(.01)/(-b[FBPFuel("C2")]);                        /* 41b - 2009 */

      if ((1. - pow(RSF_D1/a[FBPFuel("D1")], 1./c[FBPFuel("D1")])) >= 0.01)
         ISF_D1 = log(1. - pow(RSF_D1/a[FBPFuel("D1")], 1./c[FBPFuel("D1")]))
               /(-b[FBPFuel("D1")]);                                   /* 41a - 2009 */
      else
         ISF_D1 = log(.01)/(-b[FBPFuel("D1")]);                        /* 41b - 2009 */


      ISF = PC/100.*ISF_C2+(1.-PC/100.)*ISF_D1;                        /* 42a - 2009 */
   }

   if (!strcasecmp(Fuels[Fuel], "M3"))
   {
      PDF100=100.;
      RSZ = ROScalc(FBPFuel("M3"), ISZ, NoBUI, FMC, SFC, PC, PDF100, C, CBH);
      RSF_M3 = RSZ * SF;                                                    /* 40 */
      RSZ = ROScalc(FBPFuel("D1"), ISZ, NoBUI, FMC, SFC, PC, PDF, C, CBH);
      RSF_D1 = RSZ * SF;                                                    /* 40 */

      if ((1. - pow(RSF_M3/a[FBPFuel("M3")], 1./c[FBPFuel("M3")])) >= 0.01)
         ISF_M3 = log(1. - pow(RSF_M3/a[FBPFuel("M3")], 1./c[FBPFuel("M3")]))
               /(-b[FBPFuel("M3")]);                                   /* 41a - 2009 */
      else
         ISF_M3 = log(.01)/(-b[FBPFuel("M3")]);                        /* 41b - 2009 */

      if ((1. - pow(RSF_D1/a[FBPFuel("D1")], 1./c[FBPFuel("D1")])) >= 0.01)
         ISF_D1 = log(1. - pow(RSF_D1/a[FBPFuel("D1")], 1./c[FBPFuel("D1")]))
               /(-b[FBPFuel("D1")]);                                   /* 41a - 2009 */
      else
         ISF_D1 = log(.01)/(-b[FBPFuel("D1")]);                        /* 41b - 2009 */

      ISF = PDF/100.*ISF_M3 + (1.-PDF/100.)*ISF_D1;                    /* 42b - 2009 */

   }

   if (!strcasecmp(Fuels[Fuel], "M4"))
   {
      PDF100=100.;
      RSZ = ROScalc(FBPFuel("M4"), ISZ, NoBUI, FMC, SFC, PC, PDF100, C, CBH);
      RSF_M4 = RSZ * SF;                                                    /* 40 */
      RSZ = ROScalc(FBPFuel("D1"), ISZ, NoBUI, FMC, SFC, PC, PDF, C, CBH);
      RSF_D1 = RSZ * SF;                                                    /* 40 */

      if ((1. - pow(RSF_M4/a[FBPFuel("M4")], 1./c[FBPFuel("M4")])) >= 0.01)
         ISF_M4 = log(1. - pow(RSF_M4/a[FBPFuel("M4")], 1./c[FBPFuel("M4")]))
               /(-b[FBPFuel("M4")]);                                   /* 41a - 2009 */
      else
         ISF_M4 = log(.01)/(-b[FBPFuel("M4")]);                        /* 41b - 2009 */

      if ((1. - pow(RSF_D1/a[FBPFuel("D1")], 1./c[FBPFuel("D1")])) >= 0.01)
         ISF_D1 = log(1. - pow(RSF_D1/a[FBPFuel("D1")], 1./c[FBPFuel("D1")]))
               /(-b[FBPFuel("D1")]);                                   /* 41a - 2009 */
      else
         ISF_D1 = log(.01)/(-b[FBPFuel("D1")]);                        /* 41b - 2009 */

      ISF = PDF/100.*ISF_M4 + (1.-PDF/100.)*ISF_D1;                    /* 42c - 2009 */
   }

   if (!strcasecmp(Fuels[Fuel], "O1a") ||
       !strcasecmp(Fuels[Fuel], "O1b"))
   {
      if (C < 58.8)
         CF = 0.005*(exp(0.061*C)-1.);                                 /* 35a - 2009 */
      else
         CF = 0.176 + 0.02*(C-58.8);                                   /* 35b - 2009 */

      if ((1. - pow(RSF/(CF*a[Fuel]), 1./c[Fuel])) >= 0.01)
         ISF = log(1. - pow(RSF/(CF*a[Fuel]), 1.
               /c[Fuel]))/(-b[Fuel]);                                  /* 43a - 2009 */
      else
         ISF = log(0.01)/(-b[Fuel]);                                   /* 43b - 2009 */
   }

   m = 147.2*(101.-FFMC)/(59.5+FFMC);                                  /* 46 */
   fF = 91.9*exp(-.1386*m)*(1.+pow(m,5.31)/4.93e7);                    /* 45 */

//   WSE = log(ISF/(0.208 * fF))/0.05039;                                /* 44 */

   WSE = 1./0.05039*log(ISF/(0.208 * fF));                             /* 44a , 44d- 2009 */
   if (WSE > 40)                                                       /* 44e - 2009 */
      if (ISF < (0.999*2.496*fF) )
         WSE = 28.-(1./0.0818*log(1.-ISF/(2.496*fF)));                 /* 44b - 2009 */
      else
         WSE = 112.45;                                                 /* 44c - 2009 */

   WSX = WS*sin(WAZ) + WSE*sin(SAZ);                                   /* 47 */
   WSY = WS*cos(WAZ) + WSE*cos(SAZ);                                   /* 48 */
   *pWSV = sqrt(WSX*WSX + WSY*WSY);                                    /* 49 */
   *pRAZ = acos(WSY/ (*pWSV));   /* in radians */                            /* 50 */
   if (WSX < 0.) *pRAZ = 2*PI - *pRAZ;                                 /* 51 */
}

/* BUI Effect */

double BEcalc(int Fuel, double BUI)
{
double BE;
   if (BUI > 0. && BUIo[Fuel] > 0.)
      BE = exp(50.*log(q[Fuel])*(1./BUI - 1./BUIo[Fuel]));             /* 54 */
   else
      BE = 1.;

   return (BE);
}

/* Crown Fraction Burned (CFB) calculation */

double CFBcalc(int Fuel, double FMC, double SFC, double ROS, double CBH)
{
double CFB, CSI, RSO;

/* I have no idea why BMW has this in there 
   if (!strcasecmp(Fuels[Fuel], "C6")) SFC = 2.0;*/

   CFB = 0.;
   CSI = 0.001 * pow(CBH, 1.5) * pow(460. + 25.9*FMC, 1.5);            /* 56 */
   RSO = CSI/(300.*SFC);                                               /* 57 */
   if (ROS > RSO) CFB = 1. - exp(-0.23*(ROS - RSO));                   /* 58 */
   
   return (CFB);
}

void C6calc(int Fuel, double ISI, double BUI, double FMC, double SFC, double CBH, double *ROS, double *CFB, double *RSC)
{
double T, h, FME, RSI, RSS, FMEavg;

   FMEavg = 0.778;                                                     /* page 37 */
   T = 1500. - 2.75 * FMC;                                             /* 59 */
   h = 460. + 25.9 * FMC;                                              /* 60 */
   FME = pow(1.5 - 0.00275 * FMC, 4.)/(460. + 25.9*FMC) * 1000.;       /* 61 */
   RSI = 30. * pow(1. - exp(-0.08 * ISI), 3.0);                        /* 62 */
   RSS = RSI * BEcalc(Fuel, BUI);                                      /* 63 */
   *RSC = 60. * (1. - exp(-0.0497*ISI)) * FME/FMEavg;                   /* 64 */

   if (*RSC > RSS)
   {
      *CFB = CFBcalc(Fuel, FMC, SFC, RSS, CBH);
      *ROS = RSS + (*CFB)*(*RSC-RSS);                                   /* 65 */
   }
   else
   {
      *CFB=0.;
      *ROS = RSS;
   }
}

double TFCcalc(int Fuel, double CFL, double CFB, double SFC, double PC, double PDF)
{
double CFC, TFC;
   CFC = CFL * CFB;                                                    /* 66 ; 66a - 2009 */
   if (!strcasecmp(Fuels[Fuel], "M1") ||
       !strcasecmp(Fuels[Fuel], "M2"))
      CFC = PC/100.*CFC;                                               /* 66b - 2009 */
   if (!strcasecmp(Fuels[Fuel], "M3") ||
       !strcasecmp(Fuels[Fuel], "M4"))
      CFC = PDF/100.*CFC;                                              /* 66c - 2009 */
   TFC = SFC + CFC;                                                    /* 67 */

   return(TFC);
}

double FIcalc(double FC, double ROS)
{
double FI;
   FI = 300. * FC * ROS;                                               /* 69 */

   return(FI);
}

double BROScalc (int Fuel, double FFMC, double BUI, double WSV, double FMC, double SFC, double PC, double PDF, double C, double CBH)
{
double m, fF, BfW, BISI, BROS;

   m = 147.2*(101.-FFMC)/(59.5+FFMC);                                  /* 46 */
   fF = 91.9*exp(-0.1386*m)*(1.+pow(m,5.31)/4.93e7);                   /* 45 */
   BfW = exp(-0.05039*WSV);                                            /* 75 */
   BISI = 0.208*BfW*fF;                                                /* 76 */

/* Note the BUI effect is captured in ROScalc */
   BROS =  ROScalc(Fuel, BISI, BUI, FMC, SFC, PC, PDF, C, CBH);                      /* 77 */

   return (BROS);
}

double ROStcalc(int Fuel, double ROSeq, double t, double CFB)
{
double ROSt, alpha;

   if (!strcasecmp(Fuels[Fuel], "C1") ||
       !strcasecmp(Fuels[Fuel], "O1a") ||
       !strcasecmp(Fuels[Fuel], "O1b") ||
       !strcasecmp(Fuels[Fuel], "S1") ||
       !strcasecmp(Fuels[Fuel], "S2") ||
       !strcasecmp(Fuels[Fuel], "S3"))
      alpha = 0.115;                                                   /* page 41 */
   else
      alpha = 0.115 - 18.8 * pow(CFB, 2.5) * exp(-8.* CFB);            /* 72 */

   ROSt = ROSeq * (1. - exp (-alpha * t));                             /* 70 */

   return(ROSt);
}

double LBcalc(int Fuel, double WSV)
{
double LB;

   if (!strcasecmp(Fuels[Fuel], "O1a") ||
       !strcasecmp(Fuels[Fuel], "O1b"))

      if (WSV >= 1.0)
         LB = 1.1 * pow(WSV, 0.464);                                   /* corrected from "+" to "*" in the errata; 80 */
      else
         LB = 1.0;                                                     /* 81 */

   else
      LB = 1.0 + 8.729 * pow(1.-exp(-0.030*WSV), 2.155);               /* 79 */

   return(LB);
}

double LBtcalc(int Fuel, double LB, double t, double CFB) /* this is a new subroutine to account for the discussion in 3.6.1 */
{
double LBt, alpha;

   if (!strcasecmp(Fuels[Fuel], "C1") ||
       !strcasecmp(Fuels[Fuel], "O1a") ||
       !strcasecmp(Fuels[Fuel], "O1b") ||
       !strcasecmp(Fuels[Fuel], "S1") ||
       !strcasecmp(Fuels[Fuel], "S2") ||
       !strcasecmp(Fuels[Fuel], "S3"))
      alpha = 0.115;                                                   /* page 41 */
   else
      alpha = 0.115 - 18.8 * pow(CFB, 2.5) * exp(-8.* CFB);               /* 72 */

    LBt = (LB -1.) * (1. - exp(-alpha*t)) + 1.;                            /* 81 - 2009 */

   return(LBt);
}

double FROScalc(double ROS, double BROS, double LB)
{
double FROS;

   FROS = (ROS + BROS)/LB/2.;                                          /* 89 */

   return(FROS);
}

double ROSthetacalc(double ROS, double FROS, double BROS, double theta)
{
double ROStheta, c, s;

    c = cos(theta);
    s = sin(theta);
    if (c==0.)
        c = cos(theta+.001);
    else
        ROStheta = (ROS - BROS)/(2*c) + (ROS + BROS)/(2*c)*
            (FROS*c*sqrt(FROS*FROS * c*c + (ROS*BROS)*s*s) - (ROS*ROS - BROS*BROS)*s*s)/
        (    FROS*FROS*c*c + ((ROS+BROS)/2)*((ROS+BROS)/2)*s*s);           /* 94 - 2009 */

    return(ROStheta);
}

/* int Date2Julian(int month, int day, int year)
{
/* remember Jan=1, Feb=2, ... 
int days[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
int julian;

   julian = days[month-1]+day;
   if (year%4 == 0  && year > 0 && month > 2) julian = julian+1;

   return(julian); 
} */

int FBPdefaults(FBP *fbp)
{
    int status;

    strcpy(fbp->FuelType,"C2"); 
    fbp->Accel = 1;
    fbp->Dj = 0;
    fbp->Do = 0;
    fbp->ELV = -9;
    fbp->BUIEff = 0;
    fbp->t = 0;
    fbp->FFMC = 90;
    fbp->ISI = 0;
    fbp->BUI = 60;
    fbp->WS = 10;
    fbp->WD = 0;
    fbp->GS = 0;
    fbp->Aspect = 0;
    fbp->PC = 50;
    fbp->PDF = 35;
    fbp->C = 80.;
    fbp->GFL = 3.;
    fbp->CBH = 7.;
    fbp->CFL = -1.;
    fbp->LAT = 60.;
    fbp->LON = -120.;
    fbp->FMC = 0.;
    fbp->theta = 0.;
    fbp->SH = -999.;
    fbp->SD = -999.;

    status = 0;    // success
   return(status);
}

void FBPtest()
{
double SAZ, WAZ;
int TestCase, status, Date1, Date2, Do, month, day, year;

FBP *fbp;

    printf("...conducting FBP test\n");

    fbp = (FBP *)malloc(sizeof(FBP));

   infile  = fopen("fbp2009.inp", "r");
//   outfile5 = fopen("Table5.out", "w");
//   outfile6 = fopen("Table6.out", "w");
TestCase=0;

    status = FBPdefaults(fbp);

   printf("Cs Fuel    HROSt    FCFB     FROS       FFI      FROSt    BCFB     BROS       BFI      BROSt\n");

   while (!feof(infile) && TestCase<20)
   {
      fscanf(infile, "%d%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%lf%lf%lf%lf%lf", 
             &TestCase, fbp->FuelType, &fbp->FFMC, &fbp->BUI, &fbp->WS, &fbp->WD, &WAZ, 
             &fbp->GS, &SAZ, &fbp->LAT, &fbp->LON, &fbp->ELV, &fbp->Dj, &fbp->Do,
             &fbp->t, &fbp->PC, &fbp->PDF, &fbp->GFL, &fbp->C);
//     1 C1 90 130 20 0 180 15 90 55 110 -1 182 -1 20 -1 -1 -1 -1

/*      printf("%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf %lf\n\n",  
             TestCase, fbp->FuelType, fbp->FFMC, fbp->BUI, fbp->WS, fbp->WD, WAZ, 
             fbp->GS, SAZ, fbp->LAT, fbp->LON, fbp->ELV, fbp->Dj, fbp->Do,
             fbp->t, fbp->PC, fbp->PDF, fbp->GFL, fbp->C); */

      fbp->Accel = 0;
      fbp->LON = -fbp->LON;
      fbp->Aspect = SAZ-180.;
      fbp->BUIEff = 0;
      fbp->CBH = -1.;
      fbp->CFL = -1.;
      fbp->t=(fbp->t)/60.;    /* program works in decimal hours */


//      printf("%2d ", TestCase);

      status = FBPCalc(fbp);
      printf("%2d %4s %8.3lf %8.3lf %8.3lf %8.3lf %10.2lf\n", TestCase, fbp-> FuelType, fbp->SFC, fbp->CFB, fbp->TFC, fbp->ROS, fbp->HFI);
//      printf("%2d %4s %8.3lf %8.3lf %8.3lf %10.2lf %8.3lf %8.3lf %8.3lf %10.2lf %8.3lf\n", 
//          TestCase, fbp-> FuelType, fbp->HROSt,   
//          fbp->FCFB, fbp->FROS, fbp->FFI, fbp->FROSt,
//          fbp->BCFB, fbp->BROS, fbp->BFI, fbp->BROSt);

//      printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", TestCase, fbp->SFC, fbp->CFB, fbp->TFC, fbp->ROS, fbp->HFI);
//      fprintf(outfile5, "%d %lf %lf %lf %lf\n", TestCase, fbp->SFC, fbp->CFB, fbp->ROS, fbp->HFI);

   }
fclose(infile);
//fclose(outfile5);
//fclose(outfile6);

}

/* This main program is designed to read in a file and display the output to the screen 

This program can be compiled using the Gnu C compiler using the following command:

    gcc fbp_2009.c -o fbp_2009.exe

The program is run from the command prompt by passing an input file name in the command line

    fbp_2009 infile.dat

where the input file contains all the inputs in an ASCII text fomat with the variable name from the publication followed by the value, e.g.:

    FUEL C2
    FFMC 95
    ISI 15
    WS 10
    WD 180
    ...

A number of default values are used.  These are listed in lines 140-160.

Output is displayed to the screen unless the redirect command is used

    fbp_2009 infile.dat > outfile.dat
    
*/

//  To turn off this section, simply place an open comment "/*" marker prior to int main(...) and a close comment "*/" on the last line.
/*
int main(int argc, char *argv[])

{
   char infilename[64], buffer[64], errortext[100];
   double x;
   int status;
   FBP *fbp;
   
      fbp = malloc(sizeof(FBP));

//   printf("argc = %d\n", argc);

   if (argc > 1)
   {
      sprintf(infilename, "%s\0", argv[1]);
      strcpy(infilename, argv[1]);
   }
   else
   {
      printf("To run, pass an input filename to the program through the command line\n");
      printf("Using infile.dat\n");
      strcpy(infilename, "infile.dat\0");
   }
   
//    FBPtest();

    status = FBPdefaults(fbp);

    fbp->Accel = 0;


//      fbp->t=(fbp->t)/60.;    // program works in decimal hours
      
    fbp->CBH = -1.;     // NB: This was left out resulting in a default CBH=7m being used [KRA April 9, 2014]
    
//    printf("'%s'\n", infilename);

    if ( (infile = fopen(infilename, "r")) != NULL)
    {
      fscanf(infile,"%s", buffer );

      while( !feof( infile ) )
      {

         if ( !strncasecmp( buffer,"FuelType",3 ) ) 
         {
            fscanf( infile," %s", buffer);
            sprintf(fbp->FuelType, "%s\0", buffer);
         }

         if ( !strncasecmp( buffer,"Accel",3 ) ) 
            fscanf( infile," %d", &fbp->Accel);

         if ( !strcasecmp( buffer,"Dj") ) 
            fscanf( infile," %d", &fbp->Dj);

         if ( !strcasecmp( buffer,"Do" ) ) 
            fscanf( infile," %d", &fbp->Do);

         if ( !strncasecmp( buffer,"ELV",3 ) || 
              !strncasecmp( buffer, "ELEV", 3) ) 
            fscanf( infile," %d", &fbp->ELV);


         if ( !strcasecmp( buffer,"BE" ) || 
              !strncasecmp( buffer,"BUIEff" , 4 )) 
         {
            fscanf(infile," %s", buffer);
            if (!strcasecmp(buffer, "on") || !strcasecmp(buffer, "TRUE") || !strcasecmp(buffer, "yes") || atof(buffer) > 0)
                fbp->BUIEff = 1;
            else if (!strcasecmp(buffer, "off") || !strcasecmp(buffer, "FALSE") || atof(buffer) <= 0)
                fbp->BUIEff = 0;
            else
                fscanf( infile," %d", &fbp->BUIEff);
         }

         if ( !strcasecmp( buffer,"t" ) || 
              !strncasecmp( buffer, "time", 3) ) 
            fscanf( infile," %lf", &fbp->t);    // decimal hours

         if ( !strcasecmp( buffer,"FFMC") ) 
            fscanf( infile," %lf", &fbp->FFMC);

         if ( !strcasecmp( buffer,"ISI" ) ) 
            fscanf( infile," %lf", &fbp->ISI);

         if ( !strcasecmp( buffer,"BUI" ) ) 
            fscanf( infile," %lf", &fbp->BUI);

         if ( !strcasecmp( buffer,"WS" ) ) 
            fscanf( infile," %lf", &fbp->WS);

         if ( !strcasecmp( buffer,"WD" ) ||
              !strncasecmp( buffer, "WDIR",3 ) ) 
            fscanf( infile," %lf", &fbp->WD);

         if ( !strcasecmp( buffer,"WAZ" ) ) 
         {
            fscanf( infile," %lf", &fbp->WD);
            fbp->WD = fbp->WD - 180.;
         }

         if ( !strcasecmp( buffer,"GS" ) ) 
            fscanf( infile," %lf", &fbp->GS);

         if ( !strcasecmp( buffer,"Aspect" ) ) 
            fscanf( infile," %lf", &fbp->Aspect);

         if ( !strcasecmp( buffer,"SAZ" ) ) 
         {
            fscanf( infile," %lf", &fbp->Aspect);
            fbp->Aspect = fbp->Aspect - 180.;
         }

         if ( !strcasecmp( buffer,"PC" ) ) 
            fscanf( infile," %lf", &fbp->PC);

         if ( !strcasecmp( buffer,"PDF" ) ) 
            fscanf( infile," %lf", &fbp->PDF);

         if ( !strcasecmp( buffer,"C" ) ||
              !strncasecmp( buffer, "cured",3 ) ) 
            fscanf( infile," %lf", &fbp->C);

         if ( !strcasecmp( buffer,"GFL" ) ) 
            fscanf( infile," %lf", &fbp->GFL);

         if ( !strcasecmp( buffer,"CBH" ) ) 
            fscanf( infile," %lf", &fbp->CBH);

         if ( !strcasecmp( buffer,"CFL" ) ) 
            fscanf( infile," %lf", &fbp->CFL);

         if ( !strncasecmp( buffer,"LAT",3 )||
              !strcasecmp( buffer, "y" )) 
            fscanf( infile," %lf", &fbp->LAT);

         if ( !strncasecmp( buffer,"LON",3 ) ||
              !strcasecmp( buffer,"lng" ) ||
              !strcasecmp( buffer,"x" )) 
            fscanf( infile," %lf", &fbp->LON);
            

         if ( !strcasecmp( buffer,"FMC" ) ||
              !strncasecmp( buffer, "cured",3 ) ) 
            fscanf( infile," %lf", &fbp->FMC);

         if ( !strcasecmp( buffer,"SH" ) ) 
            fscanf( infile," %lf", &fbp->SH);

         if ( !strcasecmp( buffer,"SD" ) ) 
            fscanf( infile," %lf", &fbp->SD);

         if ( !strcasecmp( buffer,"theta" ) ) 
            fscanf( infile," %lf", &fbp->theta); 

//         printf("%s\n", buffer); 

         fscanf(infile,"%s", buffer );

      }
// printf("... %lf\n", &fbp->FFMC);
      fclose(infile);
      
//      fbp->t = fbp->t /60.;   // doing this only for the test data set
//      fbp->LON = -1 * fbp->LON;  // doing this only for the test data set


   status = FBPCalc(fbp);

// outputs
   printf("ROS %lf\n", fbp->ROS);              // Rate of Spread [m/min]
   printf("FROS %lf\n", fbp->FROS);            // Flank rate of Spread [m/min]
   printf("BROS %lf\n", fbp->BROS);            // Back Rate of Spread [m/min]
   if (fbp->theta > 0) printf("TROS %lf\n", fbp->TROS);            // Rate of Spread at angle theta [m/min]

   if (fbp->t > 0)
   {
      printf("HROSt %lf\n", fbp->HROSt);       // Head Rate of Spread at time t [m/min]
      printf("FROSt %lf\n", fbp->FROSt);       // Flank Rate of Spread at time t [m/min]
      printf("BROSt %lf\n", fbp->BROSt);       // Back Rate of Spread at time t [m/min]
      if (fbp->theta > 0) printf("TROSt %lf\n", fbp->TROSt);           //* Rate of Spread at angle theta at time t [m/min]
   }
   printf("CFB %lf\n", fbp->CFB);              // Crown Fraction Burned
   printf("FCFB %lf\n", fbp->FCFB);            // Flank Crown Fraction Burned [%]
   printf("BCFB %lf\n", fbp->BCFB);            // Back Crown Fraction Burned [%]
   if (fbp->theta > 0) printf("TCFB %lf\n", fbp->TCFB);            //* Crown Fraction Burned at angle thetea [%]

   printf("HFI %lf\n", fbp->HFI);              // Head Fire Intensity [kW/m]
   printf("FFI %lf\n", fbp->FFI);              // Head Fire Intensity [kW/m]
   printf("BFI %lf\n", fbp->BFI);              // Head Fire Intensity [kW/m]
   if (fbp->theta > 0) printf("TFI %lf\n", fbp->TFI);             //* Head Fire Intensity [kW/m]

   printf("TFC %lf\n", fbp->TFC);              // Total Fuel Consumption [kg/m^2]
   printf("FTFC %lf\n", fbp->FTFC);            // Flank Total Fuel Consumption [kg/m^2]
   printf("BTFC %lf\n", fbp->BTFC);            // Back Total Fuel Consumption [kg/m^2]
   if (fbp->theta > 0) printf("TTFC %lf\n", fbp->TTFC);            //* Total Fuel Consumption at angle theta [kg/m^2]

   printf("SFC %lf\n", fbp->SFC);              // Surface Fuel Consumption [kg/m^2]
   printf("TI %lf\n", fbp->TI);                // Time of Crown Fire initiation [hrs since ignition]
   printf("FTI %lf\n", fbp->FTI);              // Time of Flank Crown Fire initiation [hrs since ignition]
   printf("BTI %lf\n", fbp->BTI);              // Time of Back Crown Fire initiation [hrs since ignition]
   if (fbp->theta > 0) printf("TTI %lf\n", fbp->TTI);             //* Time of Crown Fire initiation at angle theta [hrs since ignition]

   printf("LB %lf\n", fbp->LB);                // Length to Breadth ratio
   printf("RAZ %lf\n", fbp->RAZ);              // Spread direction azimuth
   printf("WSV %lf\n", fbp->WSV);              // Net vectored wind speed

   }
   else
   {
      printf("Error reading file %s\n", infilename);
   }

   
   printf ("Done!\n\n");
return 0;
} */

int fbpcalc_(FBP *fbp)
{
    int *status;
    *status = FBPCalc(fbp);
    return(*status);
}

void fbpdefaults_(FBP *fbp)
{
    int *status;
    FBPdefaults(fbp); 
//    return(*status);
}
