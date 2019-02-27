#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "pfas.h"
#include "fbp_2009.h"
#include "feps_2011.h"

FEPS feps;
FBP fbp;
PFAS pfas;
FILE *infile, *infile2;


int ReadFEPSData(char *infilename, FEPS *feps, FBP *fbp)
{
char buffer[NB_CHAR_PATH], errortext[100];
double ndays, bias, pi = 3.1415926;
int i, status;

    status = 0;    // default = success

/* for now we are using the FBP structure and default values for CanFEPS */
    strcpy(fbp->FuelType,"C2"); 
    fbp->Accel = 1;     /* 0 = point, 1 = line (no acceleration by default unless specified; line always no acceleration) */
    fbp->Dj = 0;
    fbp->Do = 0;
    fbp->ELV = -9;
    fbp->BUIEff = -1;    /* turning off the BUI effect */
    fbp->t = 0;
    fbp->FFMC = 90;    /* This is the noon FFMC */
    fbp->ISI = 0;
    fbp->BUI = 60;
    fbp->WS = 10;
    fbp->WD = 0;
    fbp->GS = 0;
    fbp->Aspect = 0;
    fbp->PC = 50;
    fbp->PDF = 35;
    fbp->C = 80.;
    fbp->GFL = .35;
    fbp->CBH = 7.;
    fbp->CFL = -1.;
    fbp->LAT = 60.;
    fbp->LON = -120.;
    fbp->FMC = 0.;
    fbp->SH = 0.;
    fbp->SD = 0.;
    fbp->theta = 0.;
    
    fbp->CFB = 0.;
    fbp->HFI = 0.;
    fbp->TFC = 0.;
    fbp->SFC = 0.;
    
/* These are the default CFFEPS values */    

/* these are not relevant to CMC poject */
strcpy(feps->profilename, ""); /* name of Rosie's profile file name */
strcpy(feps->previous, ""); 
strcpy(feps->previous, ""); 
feps->header=-1;                /* header = off  */
feps->obs=-999;             /* target observation time (time when corresponding plume height is observed)]*/
feps->sinks=1;              /* include heat sinks in calculation >0 = yes */
feps->residencetime=0.;     /* residence time [hours] -- used to calculate the plume base size (area burned over the residence time); <=0 uses daily area */
feps->residencetime=0.;     /* residence time [hours] -- used to calculate the plume base size (area burned over the residence time); <=0 uses daily area */
feps->perimeter=-999.;      /* perimeter at dtime [m]*/


/* these are relevant to CMC poject */
    strcpy(feps->method, "cmc"); /* FIREWORK calculations */
    strcpy(feps->csvfilename, "input.csv"); /* name of Peter's hotspot file name */
    strcpy(feps->emissionsfilename, "emissions.csv"); /* name of emissions factor file name */
    strcpy(feps->emissionsmodelsfilename, "emodel.csv"); /* name of emissions model file name */
    strcpy(feps->outfilename, "output.csv"); /* name of Peter's hotspot file name */

    strcpy(feps->shape, "weighted");
    strcpy(feps->type, "average");    /* type of upper air profile method used (average lapse rate, dry, wet) */

    feps->ddate=20110701;       /* detection date of fire */
    feps->dtime=1200;           /* detection time (time all values are assumed to be collected)]*/

    feps->timezone=-999;        /* time zone - hour offset from UTC [standard time]; anything outside of -24 to 24 means timezone is calculate on the fly from the Longitude/15. */
    feps->LDT=1;                /* Local Daylight Time 1= one hour time offset */
    feps->diurnal=0;            /* include diurnal trend 0=Yes, else No */
    feps->radiation=0;          /* method of calculating radiation term: < 0 Byram's 1200/8600; =0 A_wall/A_top; > 0  value as % of Qfire */
    feps->print=0;              /* printout the data 0=FBP; 1=Q*/
    
    feps->temp=25;                /* surface temperature */
    feps->rh=40;                  /* humidity */
    feps->DMC = 50.;            /* Duff Moisture Code */
    feps->DC = 200.;            /* Drought Code */
    feps->pressure=1013.25;     /* MSL pressure [mb] */
    feps->growth=-999.;         /* area growth since last calculation [ha]*/
    feps->area=1.0;             /* size of fire at dtime [ha]*/
    feps->reset=-1;                /* reset time of smoke plume calculations [decimal hours LST] (off if < 0) */
    feps->thstart=9.;            /* time to start (exclusive) top-hat fire growth [decimal hours LST] */
    feps->thend=21.;            /* time to end (inclusive) top-hat fire growth [decimal hours LST] */

    feps->lapse=-999.;          /* lapse rate [oC/m] */
    
/* Default values taken from US standard atmosphere 1976 (http://www.digitaldutch.com/atmoscalc/) */
    feps->ts =15.;              /* temperature at surface */
    feps->t850=5.5295;             /* temperature at 850 mb */
    feps->t700=-4.578;           /* temperature at 700 mb */
    feps->t500=-21.231;          /* temperature at 500 mb */
    feps->t250=-52.3595;          /* temperature at 250 mb */
    feps->zs =111.;             /* height at surface */
    feps->z850=1457;            /* height at 850 mb */
    feps->z700=3012;            /* height at 700 mb */
    feps->z500=5574;            /* height at 500 mb */
    feps->z250=10363;            /* height at 250 mb */

    feps->alpha=12.;            /* entrainment half-angle (0. = no entrainment) */
    feps->Qo=0.;                /* amount of energy previously injected into atmosphere */
    feps->timestep=1.00;        /* timestep for modelling black body radiation loss (no timestep means no BB radiation) */
    feps->elapsed=0.;            /* elapsed time decimal hours */
    
    feps->Qfire=-999.;          /* energy of the fire: used when directly calculating the plume height from energy (in stand alone module) */
    feps->Qplume=-999.;         /* energy in the plume: used when directly calculating the plume height from energy (in stand alone module) */
    feps->Qf2Qp=-999.;          /* conversion factor from fire to plume: used when directly calculating the plume height from energy (in stand alone module) */

    feps->dz=-999.;                /* plume height [m] */
    feps->dzmin=-999.;              /* minimum plume height [m] */
            
    feps->Flaming=-999.;
    feps->Smoldering=-999.;
    feps->Residual=-999.;

   if ( (infile = fopen(infilename, "r")) != NULL)
   {
      fscanf(infile,"%s", buffer );
      while( !feof( infile ) )
      {
      
/* FBP data */
         if ( !strncasecmp( buffer,"FuelType", 3 ) ) 
         {
            fscanf( infile," %s", buffer);
            strcpy(fbp->FuelType, buffer);
          }

         if ( !strncasecmp( buffer,"Accel", 3 ) ) 
            fscanf( infile," %d", &fbp->Accel);

         if ( !strcasecmp( buffer,"Dj") ) 
            fscanf( infile," %d", &fbp->Dj);

         if ( !strcasecmp( buffer,"Do" ) ) 
            fscanf( infile," %d", &fbp->Do);

         if ( !strncasecmp( buffer,"ELV", 3 ) ||
              !strncasecmp( buffer,"ELEV", 3) ) 
            fscanf( infile," %d", &fbp->ELV);


         if ( !strcasecmp( buffer,"BE" ) || 
              !strncasecmp( buffer,"BUIEff" , 4 )) 
         {
            fscanf(infile," %s", buffer);
            if (!strcasecmp(buffer, "on") || !strcasecmp(buffer, "TRUE") || !strcasecmp(buffer, "yes") || atoi(buffer) > 0)
                fbp->BUIEff = 0;
            else if (!strcasecmp(buffer, "off") || !strcasecmp(buffer, "FALSE") || atof(buffer) <= 0)
                fbp->BUIEff = -1;
            else
                fscanf( infile," %d", &fbp->BUIEff);
         }

         if ( !strcasecmp( buffer,"FFMC") ) 
            fscanf( infile," %lf", &fbp->FFMC);

         if ( !strcasecmp( buffer,"ISI" ) ) 
            fscanf( infile," %lf", &fbp->ISI);

         if ( !strcasecmp( buffer,"BUI" ) ) 
            fscanf( infile," %lf", &fbp->BUI);

         if ( !strcasecmp( buffer,"WS" ) ||
              !strncasecmp( buffer, "SPEED", 3) ) 
            fscanf( infile," %lf", &fbp->WS);

         if ( !strcasecmp( buffer,"WD" ) ||
              !strcasecmp ( buffer, "WDIR") ||
              !strcasecmp ( buffer, "DIR") )
            fscanf( infile," %lf", &fbp->WD);

         if ( !strcasecmp( buffer,"GS" ) ||
              !strncasecmp( buffer,"SLOPE", 3 )) 
            fscanf( infile," %lf", &fbp->GS);

         if ( !strcasecmp( buffer,"Aspect" ) ) 
            fscanf( infile," %lf", &fbp->Aspect);

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
            
         if ( !strcasecmp( buffer, "ROS") ||
              !strcasecmp( buffer, "HROSt") )
           fscanf( infile," %lf", &fbp->HROSt);


        if ( !strcasecmp( buffer,"CFB" ) ) 
            fscanf( infile," %lf", &fbp->CFB);

         if ( !strcasecmp( buffer,"HFI" ) ) 
            fscanf( infile," %lf", &fbp->HFI);

         if ( !strcasecmp( buffer,"TFC" ) ) 
            fscanf( infile," %lf", &fbp->TFC);

         if ( !strcasecmp( buffer,"SFC" ) ) 
            fscanf( infile," %lf", &fbp->SFC);

/* FEPS data */
         if ( !strncasecmp( buffer,"Method", 3 ) ) 
         {
            fscanf( infile," %s", buffer);
            if (!strncasecmp(buffer, "Alberta", 3) ||
                !strcasecmp(buffer, "AB") )
                  strcpy(feps->method, "AB");
                  
            else if (!strncasecmp(buffer, "cmc", 3) ||
                !strcasecmp(buffer, "firework") )
                  strcpy(feps->method, "cmc");
                  
            else if (!strncasecmp(buffer, "test", 3))
                  strcpy(feps->method, "TEST");
                  
            else if (!strncasecmp(buffer, "Standard", 3) ||
                !strncasecmp(buffer, "ICAO", 3) )
                  strcpy(feps->method, "ICAO");
                  
            else
                  strcpy(feps->method, "cmc"); // now the default
          }

         if ( !strncasecmp( buffer,"Profilename", 4 )) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->profilename, buffer);
         }

         if ( !strncasecmp( buffer,"GEM", 3 )) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->profilename, buffer);
         }

         if ( !strncasecmp( buffer,"csvfilename", 4 )) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->csvfilename, buffer);
         }

         if ( !strncasecmp( buffer,"emissionsfilename", 4 )) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->emissionsfilename, buffer);
         }

         if ( !strncasecmp( buffer,"emodelsfilename", 4 )) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->emissionsmodelsfilename, buffer);
         }
         
         if ( !strncasecmp( buffer,"outputfilename", 3 )) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->outfilename, buffer);
         }
         
         if ( !strcasecmp( buffer,"Previous") ) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->previous, buffer);
         }

         if ( !strncasecmp( buffer,"shape", 3 ) ) 
         {
            fscanf( infile," %s", buffer);
            if (!strncasecmp(buffer, "line", 3) ||
                !strncasecmp(buffer, "wedge", 3) ||
                !strncasecmp(buffer, "perimeter", 4) )
                  strcpy(feps->shape, "line");
            else if (!strncasecmp(buffer, "crescent", 3) )
                  strcpy(feps->shape, "crescent");
            else if (!strncasecmp(buffer, "tophat", 3)  ||
                     !strncasecmp(buffer, "persistence", 4))
                  strcpy(feps->shape, "tophat");
            else if (!strncasecmp(buffer, "weighted", 4) )
                  strcpy(feps->shape, "weighted");
            else if (!strncasecmp(buffer, "ellipse", 3) )
                  strcpy(feps->shape, "ellipse");
            else
                  strcpy(feps->shape, "weighted");
          }

         if ( !strncasecmp( buffer,"type", 3 ) ) 
         {
            fscanf( infile," %s", buffer);
            if (!strncasecmp(buffer, "dry", 3) ||
                !strncasecmp(buffer, "piecewise", 4) )
                  strcpy(feps->type, "dry");        // piecewise integration of upper air profile heated to dry adiabat
            else if (!strncasecmp(buffer, "wet", 3) )
                  strcpy(feps->type, "wet");        // routine to be written
            else
                  strcpy(feps->type, "average");    // follow methodology used in Anderson et al. 2011, heating of average lapse rate to dry adiabat
          }

         if ( !strcasecmp( buffer,"ID") ) 
         {
            fscanf( infile," %s", buffer); 
            strcpy(feps->ID, buffer);
         }
         
         if ( !strncasecmp( buffer,"header" , 3) ) 
            fscanf(infile," %s", buffer);
            if (!strcasecmp(buffer, "off") || !strcasecmp(buffer, "FALSE") || !strcasecmp(buffer, "no") || atoi(buffer) != 0)
                feps->header = -1;
            else
                feps->header = 0;

         if ( !strncasecmp( buffer,"date", 3 ) ) 
         {
            fscanf( infile," %s", buffer); 
            if (strstr(buffer, "-")!=NULL)    /* assumes date is written as dd-mm-yyyy */
                feps->ddate=10000*atoi(&buffer[6])+100*atoi(&buffer[3])+atoi(buffer);
            else
                feps->ddate=atoi(buffer);
         }

         if ( !strcasecmp( buffer,"time") )
            fscanf( infile," %d", &feps->dtime);
         if ( !strcasecmp( buffer,"timezone") )
            fscanf( infile," %d", &feps->timezone);
         if ( !strncasecmp( buffer,"observation", 3 ))
            fscanf( infile," %d", &feps->obs);
         if ( !strcasecmp( buffer,"ldt" ) ) 
            fscanf( infile," %d", &feps->LDT);
         if ( !strncasecmp( buffer,"diurnal", 3 ) ) 
         {
            fscanf(infile," %s", buffer);
            if (!strcasecmp(buffer, "off") || !strcasecmp(buffer, "FALSE") || !strcasecmp(buffer, "no") || atoi(buffer) != 0)
                feps->diurnal = -1;
            else
                feps->diurnal = 0;
         }

         if ( !strncasecmp( buffer,"sinks", 3 ) ) 
            fscanf( infile," %d", &feps->sinks);
         if ( !strncasecmp( buffer,"radiation", 3 ) ) 
            fscanf( infile," %d", &feps->radiation);
         if ( !strncasecmp( buffer,"print", 5 ) ) 
            fscanf( infile," %d", &feps->print);

         if ( !strncasecmp( buffer,"temperature", 3 ) ) 
            fscanf( infile," %lf", &feps->temp); 
         if ( !strncasecmp( buffer,"rh", 2 ) ||
              !strncasecmp( buffer,"humidity", 3 )) 
            fscanf( infile," %lf", &feps->rh); 
         if ( !strcasecmp( buffer,"DMC" ) ) 
            fscanf( infile," %lf", &feps->DMC); 
         if ( !strcasecmp( buffer,"DC" ) ) 
            fscanf( infile," %lf", &feps->DC); 
         if ( !strcasecmp( buffer,"p" ) ||
              !strncasecmp( buffer,"pressure", 4 )) 
           fscanf( infile," %lf", &feps->pressure); 

         if ( !strncasecmp( buffer,"growth", 3 ) ) 
            fscanf( infile," %lf", &feps->growth);
         if ( !strncasecmp( buffer,"residencetime", 6 ) ) 
            fscanf( infile," %lf", &feps->residencetime);
         if ( !strncasecmp( buffer,"size", 3 ) ||
              !strncasecmp( buffer,"area", 3 )) 
            fscanf( infile," %lf", &feps->area);
         if ( !strncasecmp( buffer,"maxsize", 3 ) ||
              !strncasecmp( buffer,"maxarea", 3 )) 
            fscanf( infile," %lf", &feps->maxarea);
         if ( !strncasecmp( buffer,"perimeter", 3 ) ) 
            fscanf( infile," %lf", &feps->perimeter);

         if ( !strncasecmp( buffer,"reset", 3 ) ) 
            fscanf( infile," %lf", &feps->reset);

         if ( !strncasecmp( buffer,"thstart", 3 ) ) 
            fscanf( infile," %lf", &feps->thstart);
         if ( !strncasecmp( buffer,"thend", 3 ) ) 
            fscanf( infile," %lf", &feps->thend);

        if ( !strncasecmp( buffer,"lapse", 3 ) ) 
            fscanf( infile," %lf", &feps->lapse); 

         if ( !strncasecmp( buffer,"ts", 2 ) ) 
            fscanf( infile," %lf", &feps->ts); 
         if ( !strncasecmp( buffer,"t850", 2 ) ) 
            fscanf( infile," %lf", &feps->t850); 
         if ( !strncasecmp( buffer,"t700", 2 ) ) 
            fscanf( infile," %lf", &feps->t700); 
         if ( !strncasecmp( buffer,"t500", 2 ) ) 
            fscanf( infile," %lf", &feps->t500); 
         if ( !strncasecmp( buffer,"t250", 2 ) ) 
            fscanf( infile," %lf", &feps->t250); 

         if ( !strncasecmp( buffer,"zs", 2 ) ) 
            fscanf( infile," %lf", &feps->zs); 
         if ( !strncasecmp( buffer,"z850", 2 ) ) 
            fscanf( infile," %lf", &feps->z850); 
         if ( !strncasecmp( buffer,"z700", 2 ) ) 
            fscanf( infile," %lf", &feps->z700); 
         if ( !strncasecmp( buffer,"z500", 2) ) 
            fscanf( infile," %lf", &feps->z500); 
         if ( !strncasecmp( buffer,"z250", 2) ) 
            fscanf( infile," %lf", &feps->z250); 

         if ( !strncasecmp( buffer,"alpha", 3 ) ||
              !strncasecmp( buffer,"entrainment", 3 )) 
            fscanf( infile," %lf", &feps->alpha); 
         if ( !strncasecmp( buffer,"step", 3 ) ||
              !strcasecmp( buffer,"dt") ||
              !strcasecmp( buffer,"timestep") )
            fscanf( infile," %lf", &feps->timestep); 
         if ( !strncasecmp( buffer,"elapsed", 3 ) ) 
            fscanf( infile," %lf", &feps->elapsed); 

         if ( !strcasecmp( buffer,"Qo" ) ) 
            fscanf( infile," %lf", &feps->Qo); 
         if ( !strncasecmp( buffer,"Qfire", 3 ) ) 
            fscanf( infile," %lf", &feps->Qfire); 
         if ( !strncasecmp( buffer,"Qplume", 3 ) ) 
            fscanf( infile," %lf", &feps->Qplume); 
         if ( !strncasecmp( buffer,"Qf2Qp", 3 ) ) 
            fscanf( infile," %lf", &feps->Qf2Qp); 

         if ( !strncasecmp( buffer,"plume", 3 ) ||
              !strncasecmp( buffer,"height", 3 ) ||
              !strcasecmp( buffer,"dz" ) )   
            fscanf( infile," %lf", &feps->dz); 
            
         if ( !strncasecmp( buffer,"minimum", 3 ) ||
              !strcasecmp( buffer,"dzmin" ) )   
            fscanf( infile," %lf", &feps->dzmin); 

         if ( !strncasecmp( buffer,"Flaming", 3 ) ) 
            fscanf( infile," %lf", &feps->Flaming); 
         if ( !strncasecmp( buffer,"Smoldering", 3 ) ) 
            fscanf( infile," %lf", &feps->Smoldering); 
         if ( !strncasecmp( buffer,"Residual", 6 ) ) 
            fscanf( infile," %lf", &feps->Residual);

         fscanf(infile,"%s", buffer );
         
        } // while( !feof( infile ) )
      
        if (feps->obs<0)
            feps->obs = feps->dtime;

         if (feps->timezone<=-999)                   // if no timezone value,
            feps->timezone=-(int)(fbp->LON /15. + .5);    // make an educated guess at the time zone     

        if (feps->timezone>0 && fbp->LON<0)
            feps->timezone=-feps->timezone>0;
        
        if (fabs(feps->timezone)>24)
            feps->timezone=0.;

        if (!strncmp(feps->shape, "line", 3) ) 
        {
            fbp->Accel = 1; // false, no acceleration
            
        /* for a line fire, we are assigning the perimeter value to the diameter of the area [KRA 08/07/2014] */
            if (feps->area > 0.  && feps->perimeter <= 0.)
                feps->perimeter = 2.* sqrt(feps->area*10000./pi)/1000.;   /* perimeter in km */
            
        }
        if (!strncmp(fbp->FuelType, "O1", 3) ) strcpy(fbp->FuelType, "O1A");

        feps->ws = fbp->WS; /* common wind speed in both data structures */
        
        if (feps->thstart > 24. || feps->thend > 24.)
        {
            feps->thstart = (int)((int)feps->thstart / 100) + ((int)feps->thstart % 100) /60.;
            feps->thend = (int)((int)feps->thend / 100) + ((int)feps->thend % 100) /60.;
        }

        if (feps->reset > 24.)
        {
            feps->reset = (int)((int)feps->reset / 100) + ((int)feps->reset % 100) /60.;
        }
        
        if (feps->area > 0.  && feps->perimeter <= 0.) 
            feps->perimeter = 2.*pi*sqrt(feps->area*10000./pi)/1000.;   /* perimeter in km */

        if (feps->perimeter > 0. && feps->area <= 0.) 
            feps->area = 100 * 4*feps->perimeter*feps->perimeter/pi;    /* area in ha */   

        if (feps->ts ==-999.) feps->ts = feps->temp;
        if (feps->zs ==-999.) feps->zs = fbp->ELV;
        
        if (( feps->ts == 15 && feps->zs == 111) && (feps->temp != 25 || fbp->ELV != -9)) // this substitutes temp and ELV for ts and zs IFF the latter are still defaults
        {
            feps->ts = feps->temp;
            feps->zs = fbp->ELV;
        }
        
        if (feps->thstart>=feps->thend)
        {
            printf ("Warning: Top-hat start (%lf) >= Top-hat end (%lf)\n", feps->thstart, feps->thend);
            feps->thstart=9.;            /* time to start (exclusive) top-hat fire growth [decimal hours LST] */
            feps->thend=21.;            /* time to end (inclusive) top-hat fire growth [decimal hours LST] */
        }
        
        fclose(infile);        // Needs to be inside loop to avoid Segmentation Faults in UNIX

   }
   else
   {
        sprintf(errortext,"%s: Cannot open file\nProgram halted\n", infilename);
        status = 1;    // fail
   }

   if ( (infile2 = fopen(feps->previous, "r")) != NULL)
   {
      fscanf(infile2,"%s", buffer );

      while( !feof( infile2 ) )
      {
         if ( !strcasecmp( buffer,"Qo" ) ) 
            fscanf( infile2," %lf", &feps->Qo); 

         fscanf(infile2,"%s", buffer );
      }
   fclose(infile2);        // Needs to be inside loop to avoid Segmentation Faults in UNIX
   }

   return(status);
}

 int readfepsdata_(char *infilename, FEPS *feps, FBP *fbp, int *status)
{
    *status = ReadFEPSData(infilename, feps, fbp);
    return(*status);
} 

int readfepsdata2_(char *infilename, FEPS *feps, FBP *fbp)    // I can't get the return value to work
{
    int *status;
    *status = ReadFEPSData(infilename, feps, fbp);
    return(*status);
}
