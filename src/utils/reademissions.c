#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfas.h"
#include "feps_2011.h"
//#include "malloc.h"

FILE *infile;

int ReadEmissionsData(char *infilename, EMISSIONS *emissions)
{
char buffer[120], species[SPECIES_SIZE], *args[10], *token, *fields[10];
int i, iColumn, iModel, jModel, nModels, iSpecies, jSpecies, nSpecies, maxfields, iFlaming, iSmoldering, iResidual;

    for (i=0;i<MAX_SPECIES;i++)
    {
        emissions->species[i] = (char *)calloc(SPECIES_SIZE, sizeof(double));
        strncpy(emissions->species[i], "missing", SPECIES_SIZE);

/* these are the original FEPS values, smoldering and residual now replace with depth of burn (at 1 cm per hour) */
        emissions->residence_flaming=0.25;   /* residence time for flaming combustion (hours) */
        emissions->residence_smoldering=1.75; /* residence time for smoldering combustion (hours)  - follows flaming */
        emissions->residence_residual=4.;    /* residence time for residual combustion (hours)  - follows smoldering */
    }
    
    for (i=0;i<MAX_EMODELS; i++)
        emissions->emodel[i] = 0;    

    nSpecies=0; // this tracks each new species
    nModels=0;  // this tracks unique models

    
    maxfields=0;
    iSpecies = iFlaming = iSmoldering = iResidual = iModel = -999;  // these track the column number of fields within the csv (-999 = absent)
    
    for (i=0;i<10;i++)
    {
        args[i] = (char *)calloc(SPECIES_SIZE, sizeof(char));
        fields[i] =(char *)calloc(12, sizeof (char));    
    }

      if ( (infile = fopen(infilename, "r")) != NULL)
    {
        fgets(buffer, 120, infile);// header
        iColumn=0;

// read in the header
        token = strtok(buffer, ",");
        while (token != NULL && iColumn < 12)
        {
            strcpy(fields[iColumn], token);
            token = strtok(NULL, ",");

// keep track of the fields wrt columns (crude allowance for case)
            if (strstr(fields[iColumn], "SPE")) iSpecies = iColumn;
            if (strstr(fields[iColumn], "FLA")) iFlaming = iColumn;
            if (strstr(fields[iColumn], "SMO")) iSmoldering= iColumn;
            if (strstr(fields[iColumn], "RES")) iResidual = iColumn;
            if (strstr(fields[iColumn], "MOD")) iModel = iColumn;

            if (strstr(fields[iColumn], "Spe")) iSpecies = iColumn;
            if (strstr(fields[iColumn], "Fla")) iFlaming = iColumn;
            if (strstr(fields[iColumn], "Smo")) iSmoldering= iColumn;
            if (strstr(fields[iColumn], "Res")) iResidual = iColumn;
            if (strstr(fields[iColumn], "Mod")) iModel = iColumn;

            if (strstr(fields[iColumn], "spe")) iSpecies = iColumn;
            if (strstr(fields[iColumn], "fla")) iFlaming = iColumn;
            if (strstr(fields[iColumn], "smo")) iSmoldering= iColumn;
            if (strstr(fields[iColumn], "res")) iResidual = iColumn;
            if (strstr(fields[iColumn], "mod")) iModel = iColumn;
            
            iColumn++;
        }

        maxfields=iColumn;

// read in the rest of the file

        fgets(buffer, 120, infile);// first line of data

         while( !feof( infile ) && nModels < MAX_EMODELS && nSpecies < MAX_SPECIES)
        {
            iColumn=0; 
            
            token = strtok(buffer, ",");

            while (token != NULL&& iColumn < maxfields)
            {
                strcpy(args[iColumn], token);
                token = strtok(NULL, ",");
                iColumn++;
            }
            
// Find jModel 
            if (iModel==-999)   // there is no column for emissions model
            {
                jModel=0;
                nModels=0;
            }
            else
            {
                jModel = -1;
                for (i=0; i<nModels; i++)
                    if (emissions->emodel[i] == atoi(args[iModel])) jModel = i;
            
                if (jModel<0)  // new model
                {
                    jModel=i;
                    emissions->emodel[jModel] = atoi(args[iModel]);  // add to stack of models within the list
                    nModels++;
                }
            }

// Find species
            jSpecies = -1;
            for (i=0; i<nSpecies; i++)
                if (!strcmp(emissions->species[i], args[iSpecies])) jSpecies = i;

            if (jSpecies<0)
            {
                jSpecies=i;
                strcpy(emissions->species[jSpecies], args[iSpecies]);  // add to stack of species within the list
                nSpecies++;
            }
            
            if (iFlaming >= 0)    emissions->a[jModel][jSpecies] = atof(args[iFlaming]);         else emissions->a[jModel][jSpecies] = 0; 
            if (iSmoldering >= 0) emissions->b[jModel][jSpecies] = atof(args[iSmoldering]);      else emissions->b[jModel][jSpecies] = 0; 
            if (iResidual >= 0)   emissions->c[jModel][jSpecies] = atof(args[iResidual]);        else emissions->c[jModel][jSpecies] = 0;

// printf("... %d %d %d %s %lf %lf %lf\n", jModel, jSpecies, emissions->emodel[jModel], emissions->species[jSpecies], emissions->a[jModel][jSpecies] ,emissions->b[jModel][jSpecies], emissions->c[jModel][jSpecies]);

            fgets(buffer, 120, infile);
        }
        fclose(infile);

        emissions->nModels=nModels;
        emissions->nSpecies=nSpecies;
        
        return (nSpecies);
    }
    else
        return(-1);
}


int ReadEmissionsModels(char *infilename, EMODELS *emodels)
{
    char buffer[120];
    char *token;
    int i, j, k, found;
    
    for (i=0;i<MAX_CWFIS;i++)
        emodels->FuelType[i] = (char *)calloc(4, sizeof(char));
    
// printf("Reading emissions model data: %s...\n", infilename);

    j=0; emodels->nModels = 0;
    
      if ( (infile = fopen(infilename, "r")) != NULL)
    {
        fgets(buffer, 120, infile);// header

        fgets(buffer, 120, infile);// first line of data

        while( !feof( infile ) && j < MAX_CWFIS)
        {
            token = strtok(buffer, ",");
            strncpy(emodels->FuelType[j], token, 3);
            token=strtok(NULL, ",");
            emodels->emodel[j] = atoi(token);
            
            found = 0;
            for (i=0; i<j; i++)
                if (emodels->emodel[i] == emodels->emodel[j]) found = 1;
            
            if (found==0) emodels->nModels++;

// printf("%s %d %d\n", emodels->FuelType[j], emodels->emodel[j], emodels->nModels);
            fgets(buffer, 120, infile);
            j++;
        }
        return(j);
    }
    else
        return(-1);
}


int reademissionsdata_(char *infilename, EMISSIONS *emissions, int *status)
{
//    int *status;
    *status = ReadEmissionsData(infilename, emissions);
    return(*status);
}


int reademissionsmodels_(char *infilename, EMODELS *emodels, int *status)
{
//    int *status;
    *status = ReadEmissionsModels(infilename, emodels);
    return(*status);
}

