#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
/*#include "fwi84.h" */

FILE *infile;
     
/* FFMC calculation **************************************************/

double FFMCcalc(T,H,W,ro,Fo)
double T,
       H,
       W,
       ro,
       Fo;

{
  double mo,rf,mr,Ed,Ew,ko,kd,kl,kw,m,F;
          
     mo = 147.2*(101.-Fo)/(59.5+Fo);                       /* 1  */
     if(ro > 0.5) 
        {
        rf = ro - 0.5;                                      /* 2  */
        if(mo <= 150.)
           mr = mo + 
           42.5*rf*(exp(-100./(251.-mo)))*(1-exp(-6.93/rf)); /* 3a */
        else
           mr = mo + 
                42.5*rf*(exp(-100./(251.-mo)))*(1-exp(-6.93/rf)) +
                .0015*pow(mo-150.,2.)*pow(rf,.5);            /* 3b */
        if(mr > 250.)
           mr = 250.;
        mo = mr;
        }
     Ed = 0.942*pow(H,.679)
          + 11.*exp((H-100.)/10.)+.18*(21.1-T)*(1.-exp(-.115*H));   /* 4  */
     if(mo > Ed)
        {
        ko = 0.424*(1.-pow(H/100.,1.7))
             + 0.0694*pow(W,.5)*(1.-pow(H/100.,8.));       /* 6a */
        kd = ko*.581*exp(0.0365*T);                       /* 6b */
        m = Ed+(mo-Ed)*pow(10.,-kd);                         /* 8  */
        }
     else
        {
        Ew = 0.618*pow(H,.753)
          + 10.*exp((H-100.)/10.)
          + .18*(21.1-T)*(1.-exp(-.115*H));              /* 5  */
        if(mo < Ew)
           {
           kl = 0.424*(1.-pow((100.-H)/100.,1.7))
            + 0.0694*pow(W,.5)*(1-pow((100.-H)/100.,8.));  /* 7a */
           kw = kl*.581*exp(0.0365*T);                    /* 7b */
           m = Ew - (Ew - mo)*pow(10.,-kw);                  /* 9  */
           }
        else
           m = mo;
        }
     F = 59.5*(250.-m)/(147.2+m);                            /* 10 */

  return (F);
}

/* DMC calculation ***************************************************/

double DMCcalc(T,H,ro,Po,I)
double T,
       H,
       ro,
       Po;
int    I;

{     
  double re,Mo,Mr,K,b,P,Pr;
  static double Le[] = {6.5,7.5,9.,12.8,13.9,13.9,12.4,10.9,9.4,8.,7.,6.};
     
     if(ro > 1.5)
        {
        re = 0.92*ro-1.27;                                    /* 11  */
        Mo = 20. + exp(5.6348 - Po/43.43);                    /* 12  */
        if(Po <= 33.)
           b = 100./(.5+.3*Po);                               /* 13a */
        else 
           if(Po <= 65.)
              b = 14. - 1.3*(log(Po));                        /* 13b */
           else
              b = 6.2*log(Po) - 17.2;                         /* 13c */
        Mr = Mo +1000.*re/(48.77+b*re);                       /* 14  */
        Pr = 244.72 - 43.43*log(Mr - 20.);                    /* 15  */
        if(Pr > 0.) 
           Po = Pr;
        else
           Po = 0.;
        }
     if(T > -1.1)
        K = 1.894*(T +1.1)*(100.-H)*Le[I-1]*1.0E-6;           /* 16  */
     else
        K = 0.;
     P = Po+100.*K;                                           /* 17  */

  return(P);
}

/* DC calculation ****************************************************/

double DCcalc(T,ro,Do,I)
double T,
       ro,
       Do;
int    I;

{
     
  double rd,Qo,Qr,V,D,Dr;
  static double Lf[] = {-1.6,-1.6,-1.6,.9,3.8,5.8,6.4,5.,2.4,.4,-1.6,-1.6};

     if(ro > 2.8)
        {
        rd = 0.83*(ro) - 1.27;                                /* 18  */
        Qo = 800.*exp(-Do/400.);                              /* 19  */
        Qr = Qo +3.937*rd;                                    /* 20  */
        Dr = 400.*log(800./Qr);                               /* 21  */
        if(Dr > 0.)
           Do = Dr;
        else
           Do = 0.;
        }
     if(T > -2.8)
        V = 0.36*(T+2.8)+Lf[I-1];                             /* 22  */
     else
        V = Lf[I-1];
     if(V < 0.)
        V = 0.;
     D = Do + 0.5*V;                                          /* 23  */

  return(D);
}

/* ISI calculation ***************************************************/

double ISIcalc(F,W)
double F,
       W;

{
  double fW,m,fF,R;
     fW = exp(0.05039*W);                                  /* 24  */
     m = 147.2*(101-F)/(59.5+F);                              /* 1   */
     fF = 91.9*exp(-.1386*m)*(1.+pow(m,5.31)/4.93E7);         /* 25  */
     R = 0.208*fW*fF;                                         /* 26  */
     
  return(R);
}

/* BUI calculation ***************************************************/

double BUIcalc(P,D)
double P,
       D;
{
  double U;

    if (P*D == 0)         // updated fix 2015-11-09 KRA 
        U=0;
    else
  
        if(P <= .4*D)
            U = 0.8*P*D/(P+.4*D);                                 /* 27a */
        else
            U = P - (1.-.8*D/(P+.4*D))
                *(.92+pow(.0114*P,1.7));                           /* 27b */
  return(U);
}

/* FWI calculation ***************************************************/

double FWIcalc(R,U)
double R,
       U;
{
     
  double fD,B,S;

     if(U <= 80.)
        fD = .626*pow(U,.809)+2.;                             /* 28a */
     else
        fD = 1000./(25.+108.64*exp(-.023*U));                 /* 28b */
     B = .1*R*fD;                                             /* 29  */
     if(B > 1.)
        S = exp(2.72*pow(.434*log(B),.647));                  /* 30a */
     else
        S = B;                                                /* 30b */

  return(S);
}

#include <math.h>

/* DSR calculation ***************************************************/

double DSRcalc(S)
double S;
{
   double DS;

      DS = .0272*pow(S,1.77);                               /* 31  */

   return(DS);
}

/********************** Kerry Anderson's Stuff ***************************/

/* This int main() section provides a means of calculating FWI values from the command line.

This program can be compiled using the Gnu C compiler using the following command:

    gcc fwi84.c -o fwi84.exe
    
The program is run from the command prompt by passing an input file name in the command line

    fwi84 infile.dat
    
Required inputs which are stored in an input file.  

The name of the input file is passed to the program through the command line.

Input values are passed, one per line, with a Key word followed by a value:

      Noon Weather Values:
       Key   Value
         T   Temperature (oC)
         H   Humidity (%)
         W   Wind Speed (km/hr)
         ro  Past 24 hour rainfall (mm)
         I   Month (Jan=1, etc.)

      Indices:
        Key  Value
         Fo  Yesterday's Fine Fuel Moisture Code
         Po  Yesterday's Duff Moisture Code
         Do  Yesterday's Drought Code

Key words are case insensitive, use either the above Key words (T, H, W, ...)  or reasonable proximieties (e.g. FFMC, Temp, RH, etc.)

A carriage return is required after the last value.

Default values are embedded within the code (starting values, best guess weather) to avoid crashes, 
but it is wise to input ALL values.

Output is displayed to the screen unless the redirect command is used

    fwi84 infile.dat > outfile.dat


*/

//  To turn off this section, simply place an open comment "/*" marker prior to int main(...) and a close comment "*/" on the last line.
/*

int main(int argc, char *argv[])
{
   char infilename[64], buffer[64], errortext[100], month[64]; 
   double x, T, H, W, ro, Fo, Po, Do, F, P, D, R, U, S, DSR;
   int i, I, status;
   static char *months[]=  {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

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

// Default values
   T = 25.0;
   H = 40.;
   W = 10.;
   ro = 0.;
   I = 7;
   Fo = 85.;
   Po = 6.;
   Do = 15.;
   
   if ( (infile = fopen(infilename, "r")) != NULL)
   {
      fscanf(infile,"%s", buffer );

      while( !feof( infile ) )
      {

         if ( !strncasecmp( buffer,"TEMP", 3) ||
              !strcasecmp( buffer,"T" ) ) 
            fscanf( infile," %lf", &T);

         if ( !strcasecmp( buffer,"RH" ) ||
              !strncasecmp( buffer,"HUMIDITY", 3) ||
              !strcasecmp( buffer,"H" ) ) 
            fscanf( infile," %lf", &H);

         if ( !strncasecmp( buffer,"WIND", 3 ) ||
              !strcasecmp( buffer,"WS" ) ||
              !strcasecmp( buffer,"W" ) ) 
            fscanf( infile," %lf", &W);

         if ( !strncasecmp( buffer,"RAIN", 3 ) ||
              !strncasecmp( buffer,"PRECIP", 3 ) ||
              !strcasecmp( buffer,"ro" ) ||
              !strcasecmp( buffer,"R" ) ) 
            fscanf( infile," %lf", &ro);

         if (!strncasecmp( buffer,"MONTH", 3) )
         {
            fscanf (infile, " %s", month);
            if (atoi(month) > 0 && atoi(month)<13)
               I= atoi(month);
            else
               for (i=0; i<12; i++)
                  if (!strncasecmp (month, months[i], 3) ) I = i+1;
           
          }  

         if ( !strcasecmp( buffer,"FFMC") ) 
            fscanf( infile," %lf", &Fo);

         if ( !strcasecmp( buffer,"DMC" ) ) 
            fscanf( infile," %lf", &Po);

         if ( !strcasecmp( buffer,"DC" ) ) 
            fscanf( infile," %lf", &Do);

         fscanf(infile,"%s", buffer );
        }
    }
    
 // Token limits on input weather
    if (T < -50. || T > 50.) T = 25.;
    if (H < 0. || H > 100.) H = 40.;
    if (W < 0. || W > 200.) W = 10.;
    if (ro < 0. || ro > 500.) ro = 0.;
    if (Fo < 0. || Fo > 101.) F = 85.;
    if (Po < 0. || Po > 300.) P = 6.;
    if (Do < 0. || Do > 1200.) D = 15.;
    

    F = FFMCcalc(T,  H,  W,  ro,  Fo);
    P = DMCcalc(T,  H,  ro,  Po,  I);
    D = DCcalc(T,  ro,  Do,  I);
    R = ISIcalc(F,  W);
    U = BUIcalc(P,  D);
    S = FWIcalc(R,  U);
    DSR = DSRcalc(S);
    
    printf("FFMC %lf\n", F);
    printf("DMC %lf\n", P);
    printf("DC %lf\n", D);
    printf("ISI %lf\n", R);
    printf("BUI %lf\n", U);
    printf("FWI %lf\n", S);
    printf("DSR %lf\n", DSR);

}
*/
