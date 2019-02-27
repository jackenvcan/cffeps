/*------------------------------------------------------------------
  Program..: FBP97PRO.C
  Author...: Kerry Anderson, ace programmer
  Date.....: 1Jan 1, 1993
  Notice...: Copyright (c) 1993, Forestry Canada, All Rights Reserved
  Notes....: FBP 91 fire shape parameter calculations
-------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fbp_2009.h"

/* function prototypes */
double converge(double, double, int);

double converge(alpha,T,Accel)
double alpha, T;
int Accel;
{
double t;
int i;
   t = T;
   if (!Accel)
      for (i=0;i < 10;i++)
        t = T - exp(-alpha*t)/alpha + 1/alpha;                         /* 71 */
   t = t/60;
   return(t);
}

void PROCalc(Fuel,Accel,CFB,Rh,Rf,Rb,t1,t2,A,P,D)
char   *Fuel;       /* The only significant fuel type is 4 (C1)        */
int    Accel;       /* Are Acceleration effects be incorporated? ( 0 = True (point), 1 = False (line))*/
double CFB,         /* Crown Fraction Burned                           */
       Rh,          /* Headfire Eq. Rate of spread in m/min            */
       Rf,          /* Flankfire Eq. Rate of spread in m/min           */
       Rb,          /* Backfire Eq. Rate of spread in m/min            */
       *t1,         /* Elapsed time 1 in hrs                           */
       *t2,         /* Elapsed time 2 in hrs                           */
       *A,          /* Area burned in heactares                        */
       *P,          /* Perimeter encompassed in kilometres             */
       *D;          /* Distance travelled in kilometres                */

/*
      From one non-zero fire shape input parameter (area, perimeter, forward 
      spread distance, time1, time2), PROCalc calculates all other fire shape
      parameters.  Accelleration effects may be included.

      Note that all desired outputs must be filled with 0 when called. If all
      input parameters are 0, the outputs are zero.

      If the headfire rate of spread (Rh) is zero, then the second elapsed time
      is set to zero and all remaining parameters are untouched (presumably
      all but one have already been set to zero by the calling routine).
      
*/
{
   double pi,alpha,T,T1,T2,a,b,M;

   pi = 3.1415926;
   if (!strcmp(Fuel, "C1") || CFB==0.)
      alpha = 0.115;
   else
      alpha = (0.115 - 18.8 * pow(CFB,2.5) * exp(-8.* CFB));           /* 72 */

   if(Rh != 0.)                                           /* rate of spread */
   {
      a = (Rh+Rb)/2;
      b = Rf;
      M = (a-b)/(a+b);
   }

   if(Rh == 0.)                                        /* no rate of spread */
   {
      *t2 = 999.9;
   }
   else if(*t2 > 0.)                                                   /* hrs */
   {
      if (!Accel) 
      {
         T2 = *t2*60;
         T1 = *t1*60;
      }
      else
      {
         T2=(*t2*60+(exp(-alpha* *t2*60)-1)/alpha);                    /* 71a */
         T1=(*t1*60+(exp(-alpha* *t1*60)-1)/alpha);                    /* 71a */
      }

/* NB: values calculated are the changes in distance/area/perimeter 
       over the time interval */

      *D = (Rh*T2)/1000
         - (Rh*T1)/1000;                                               /* 71 */
      *A = (pi/2*(Rh + Rb)* Rf*T2*T2)/10000
         - (pi/2*(Rh + Rb)* Rf*T1*T1)/10000;                           /* 84 */
      *P = (pi*T2*((Rh+Rb)/2+Rf)*(1+M*M/4))/1000
         - (pi*T1*((Rh+Rb)/2+Rf)*(1+M*M/4))/1000;                      /* 85 */
   }
   else if(*A > 0.)                                                    /* ha  */
   {
      T = sqrt(*A/(pi/2*(Rh + Rb)* Rf))*100;                           /* 84 */
      *D = (Rh*T)/1000;                                                /* 71 */
      *P = (pi*T*((Rh+Rb)/2+Rf)*(1+M*M/4))/1000;                       /* 85 */
      *t2 = converge(alpha,T,Accel);
      *t1 = 0.;
   }
   else if(*P > 0.)                                                    /* km  */
   {
      T = *P/(pi*((Rh+Rb)/2+Rf)*(1+M*M/4))*1000;                       /* 85 */
      *D = (Rh*T)/1000;                                                /* 71 */
      *A = (pi/2*(Rh + Rb)* Rf*T*T)/10000;                             /* 84 */
      *t2 = converge(alpha,T,Accel);
      *t1 = 0.;
   }
   else if(*D > 0.)                                                    /* km  */
   {
      T = (*D/ Rh)*1000;                                               /* 71 */
      *A = (pi/2*(Rh + Rb)* Rf*T*T)/10000;                             /* 84 */
      *P = (pi*T*((Rh+Rb)/2+Rf)*(1+M*M/4))/1000;                       /* 85 */
      *t2 = converge(alpha,T,Accel);
      *t1 = 0.;
   }
}
