#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pfas.h"
#include "fbp_2009.h"
#include "feps_2011.h"

FEPS feps;
FBP fbp;
PFAS pfas;
FILE *fepsoutfile;


int WriteFEPSData(char *fepsoutfilename, FEPS *feps, FBP *fbp)
{
char buffer[64], errortext[100];
double ndays, bias, pi = 3.1415926;
int status;

    status = 0;    // default = success

    if ( (fepsoutfile = fopen(fepsoutfilename, "w")) != NULL)
    {
    
    /* FBP data */
        fprintf(fepsoutfile, "FuelType %s\n", fbp->FuelType);    
        fprintf(fepsoutfile, "Accel %d\n", fbp->Accel);
        fprintf(fepsoutfile, "Dj %d\n", fbp->Dj);
        fprintf(fepsoutfile, "Do %d\n", fbp->Do);
        fprintf(fepsoutfile, "Elev %d\n", fbp->ELV);
        fprintf(fepsoutfile, "BUIEff %d\n", fbp->BUIEff);    /* turning off the BUI effect */
        fprintf(fepsoutfile, "FFMC %lf\n", fbp->FFMC);    /* This is the noon FFMC as read into FGM, but the diurnally adjusted in ENERGYx? */
        fprintf(fepsoutfile, "ISI %lf\n", fbp->ISI);
        fprintf(fepsoutfile, "BUI %lf\n", fbp->BUI);
        fprintf(fepsoutfile, "WS %lf\n", fbp->WS);
        fprintf(fepsoutfile, "WD %lf\n", fbp->WD);
        fprintf(fepsoutfile, "GS %lf\n", fbp->GS);
        fprintf(fepsoutfile, "Aspect %lf\n", fbp->Aspect);
        fprintf(fepsoutfile, "PC %lf\n", fbp->PC);
        fprintf(fepsoutfile, "PDF %lf\n", fbp->PDF);
        fprintf(fepsoutfile, "C %lf\n", fbp->C);
        fprintf(fepsoutfile, "GFL %lf\n", fbp->GFL);
        fprintf(fepsoutfile, "CBH %lf\n", fbp->CBH);
        fprintf(fepsoutfile, "CFL %lf\n", fbp->CFL);
        fprintf(fepsoutfile, "LAT %lf\n", fbp->LAT);
        fprintf(fepsoutfile, "LNG %lf\n", fbp->LON);
        fprintf(fepsoutfile, "FMC %lf\n", fbp->FMC);
        fprintf(fepsoutfile, "SH %lf\n", fbp->SH);
        fprintf(fepsoutfile, "SD %lf\n", fbp->SD);
        fprintf(fepsoutfile, "THETA %lf\n", fbp->theta);

        fprintf(fepsoutfile, "HROSt %lf\n", fbp->HROSt);
        fprintf(fepsoutfile, "CFB %lf\n", fbp->CFB);
        fprintf(fepsoutfile, "HFI %lf\n", fbp->HFI);
        fprintf(fepsoutfile, "TFC %lf\n", fbp->TFC);
        fprintf(fepsoutfile, "SFC %lf\n", fbp->SFC);

        
    /* FEPS data */
        fprintf(fepsoutfile, "Method %s\n", feps->method);
        if (strlen(feps->previous) > 0) fprintf(fepsoutfile, "Previous %s\n", feps->previous);
        fprintf(fepsoutfile, "shape %s\n", feps->shape);
        fprintf(fepsoutfile, "ID %s\n", feps->ID);
        
        fprintf(fepsoutfile, "header %d\n", feps->header);
        fprintf(fepsoutfile, "date %d\n", feps->ddate);
        fprintf(fepsoutfile, "time %d\n", feps->dtime);
        fprintf(fepsoutfile, "ldt %d\n", feps->LDT);
        fprintf(fepsoutfile, "diurnal %d\n", feps->diurnal);
        fprintf(fepsoutfile, "sinks %d\n", feps->sinks);
        fprintf(fepsoutfile, "radiation %d\n", feps->radiation);
        fprintf(fepsoutfile, "print %d\n", feps->print);
          
        fprintf(fepsoutfile, "temperature %lf\n", feps->temp);
        fprintf(fepsoutfile, "rh %lf\n", feps->rh);
//        fprintf(fepsoutfile, "ws %lf\n", feps->ws);  common to fbp
        fprintf(fepsoutfile, "DMC %lf\n", feps->DMC);
        fprintf(fepsoutfile, "DC %lf\n", feps->DC);
        fprintf(fepsoutfile, "pressure %lf\n", feps->pressure);
        
        fprintf(fepsoutfile, "growth %lf\n", feps->growth);
        fprintf(fepsoutfile, "area %lf\n", feps->area);
        fprintf(fepsoutfile, "perimeter %lf\n", feps->perimeter);
        
        fprintf(fepsoutfile, "lapse %lf\n", feps->lapse);
        fprintf(fepsoutfile, "ts %lf\n", feps->ts);
        fprintf(fepsoutfile, "t850 %lf\n", feps->t850);
        fprintf(fepsoutfile, "t700 %lf\n", feps->t700);
        fprintf(fepsoutfile, "t500 %lf\n", feps->t500);
        
        fprintf(fepsoutfile, "zs %lf\n", feps->zs);
        fprintf(fepsoutfile, "z850 %lf\n", feps->z850);
        fprintf(fepsoutfile, "z700 %lf\n", feps->z700);
        fprintf(fepsoutfile, "z500 %lf\n", feps->z500);
        
        fprintf(fepsoutfile, "alpha %lf\n", feps->alpha);
        fprintf(fepsoutfile, "timestep %lf\n", feps->timestep);
        fprintf(fepsoutfile, "elapsed %lf\n", feps->elapsed);

        fprintf(fepsoutfile, "Qo %lf\n", feps->Qo);
        fprintf(fepsoutfile, "Qfire %lf\n", feps->Qfire);
        fprintf(fepsoutfile, "Qplume %lf\n", feps->Qplume);
        fprintf(fepsoutfile, "Qf2Qp %lf\n", feps->Qf2Qp);
        
        fprintf(fepsoutfile, "dz %lf\n", feps->dz);
        
        fprintf(fepsoutfile, "Flaming %lf\n", feps->Flaming);
        fprintf(fepsoutfile, "Smoldering %lf\n", feps->Smoldering);
        fprintf(fepsoutfile, "Residual %lf\n", feps->Residual);
        
        status = 0;

        fclose (fepsoutfile);   
    }
    else
   {
        sprintf(errortext,"%s: Cannot open file\nProgram halted\n", fepsoutfilename);
        status = 1;
   }

   return(status);
    
}

int writefepsdata_(char *fepsoutfilename, FEPS *feps, FBP *fbp, int *status)
{
    *status = WriteFEPSData(fepsoutfilename, feps, fbp);
    return(*status);
} 
