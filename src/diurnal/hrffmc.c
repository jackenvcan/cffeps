/*
From mwotton@pnfi.forestry.ca  Tue Mar 23 10:43:02 1993
Date: Tue, 23 Mar 93 12:39:30 EST
From: mwotton@pnfi.forestry.ca (Mike Wotton)
To: kanderson@nofc.forestry.ca
Subject: c code four hourly ffmc stuff
Content-Length: 1095

kerry
  rob tells me you were looking for this stuff....enjoy...let me know if
you have any problems
bmw 
*/

/* modified Dec 23, 1996 for data type consistency */

#include <math.h>

double hourly_ffmc(double temp, double rh, double wind, double rain, double oldffmc)
{
  double rf=42.5,drf=0.0579;
  double mo,ed,ew,moew,moed,xm,a1,e,moe,xkd;

  if (oldffmc < 0 || oldffmc > 101 ||
      temp < -50. || temp > 50 ||
      rh < 0      || rh > 100 ||
      wind < 0    || wind > 200 ||
      rain < 0    || rain > 200)
      return (-999.);
  else
  {
     mo=147.2*(101.-oldffmc)/(59.5+oldffmc);

     if(rain!=0.)
     {
        mo+=rain*rf*exp(-100./(251.-mo))*(1.0-exp(-6.93/rain));
        if(mo>250.) mo=250.;
     }

     ed=0.942*pow(rh,0.679)+(11.0*exp( (rh-100.)/10.))+0.18*(21.1-temp)*
        (1.-1./exp(0.115*rh));

     moed=mo-ed;
     ew=0.618*pow(rh,0.753)+(10.*exp((rh-100.)/10.))+0.18*(21.1-temp)*
        (1.-1./exp(0.115*rh));

     moew=mo-ew;

     if (moed==0. || (moew>=0. && moed<0.))
     {
        xm=mo;
        if(moed==0.) e=ed;
        if(moew>=0.) e=ew;
     }
     else
     {
        if( moed>0.)
        {
           a1=rh/100.;
           e=ed;
           moe=moed;
        }
        else
        {
           a1=(100.-rh)/100.;
           e=ew;
           moe=moew;
        }

        xkd=(0.424*(1.-pow(a1,1.7))+(0.0694*sqrt(wind)*(1.-pow(a1,8.))));
//        drf = 0.581/24;
        xkd=xkd*drf*exp(0.0365*temp);
        xm=e+moe*exp(-2.303*xkd);
//        xm=e+moe*pow(10.,-xkd);
     }

     mo=xm;
     return ( 59.5*(250.-xm)/(147.2+xm) );
  }
}


