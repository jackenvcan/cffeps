/****************************************************************************/
/****************       Diurnal FFMC Adjustment Program      ****************/
/*                                                                          */
/* Date: March 2, 1993                                                      */
/*                                                                          */
/* Authors: Brad Armitage                       Bruce Lawson                */
/*          Fire Research Technician            Head, Fire Research Program */
/*          Pacific Forestry Centre             Pacific Forestry Centre     */
/*          506 W. Burnside Rd.                 506 W. Burnside Rd.         */
/*          Victoria, B.C.                      Victoria, B.C.              */
/*          Phone (604) 363-0693                Phone (604) 363-0710        */
/*          E-mail barmitage@pfc.forestry.ca                                */
/*              or   blawson@pfc.forestry.ca                                */ 
/*                                                                          */
/* Slightly modified by Kerry Anderson, programmer extraorinaire:           */
/*                                                                          */
/* 1. Calc_New_FFMC accepts a negative rh as missing (90 is substituted in) */
/* 2. routines main, InputError, Input3_Output and Input4_Output            */
/*    have been removed                                                     */
/****************************************************************************/

#include <math.h>
#include <stdio.h>

/*------------------------- Curve Equations ----------------------------------*/

/*--------------------------------------------------------------------------*/
double interp(double current, double next, double min, double timestep)
/*--------------------------------------------------------------------------*/
{ double slope = 0.0;

  slope = ((next-current)/timestep);
  return(slope*min);
}

/*--------------------------------------------------------------*/
double g31200a(double x)
/*--------------------------------------------------------------*
   Eqn# 7113  ûy=(a+cx+exý)/(1+bx+dxý)
 *--------------------------------------------------------------*/
{
  double y;
  y=(1.460075955727318+x*(0.2815668298889479+
    x*-0.01282068681081192))/
    (1.0+x*(-0.0003907916679793910+x*-0.001539833829955819));
  return(y*y);
}

/*--------------------------------------------------------------*/
double g31200b(double x)
/*--------------------------------------------------------------*
   Eqn# 4740  y=a+bx+cx3+dûx+eexp-x
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x*x;
  x3=sqrt(x);
  x4=exp(-x);
  y= -60.05817857808971-0.7922650736402593*x1
    +1.049360630941214E-05*x2+24.04228773333402*x3
    -4790555393.429649*x4;
  return(y);
}


/*--------------------------------------------------------------*/
double g31300a(double x)
/*--------------------------------------------------------------*
   Eqn# 7114  ûy=(a+cx+exý)/(1+bx+dxý+fx3)
 *--------------------------------------------------------------*/
{
  double y;
  y=(1.255216373497914+x*(0.3580951800872066+
    x*-0.01642423056280615))/
    (1.0+x*(0.02292170670740257+x*(-0.003331105559008215+
    x*3.056635174762784E-05)));
  return(y*y);
}

/*--------------------------------------------------------------*/
double g31300b(double x)
/*--------------------------------------------------------------*
   Eqn# 4560  y=a+bx+cxýlnx+d/ûx+elnx/x
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x*log(x);
  x3=1.0/sqrt(x);
  x4=log(x)/x;
  y=806.4657627391588-1.491623456269165*x1
    +0.0008873186770365859*x2-11465.74581162260*x3
    +12093.78039899465*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31400a(double x)
/*--------------------------------------------------------------*
   Eqn# 4599  y=a+bx+cxýûx+dexpx+elnx
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x*sqrt(x);
  x3=exp(x);
  x4=log(x);
  y=0.9082173869418283+0.9897247521440231*x1
    +0.001041606201107530*x2+4.633997000563253E-11*x3
    -0.005581970517308027*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31400b(double x)
/*--------------------------------------------------------------*
   Eqn# 4874  y=a+bx+cûxlnx+dx/lnx+e/xý
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=sqrt(x)*log(x);
  x3=x/log(x);
  x4=1.0/(x*x);
  y=6403.107752693009+352.7042531428901*x1
    +873.3642943916943*x2-3766.492574730166*x3
    +3580.933366117637*x4;
  return(y);
}


/*--------------------------------------------------------------*/
double g31500a(double x)
/*--------------------------------------------------------------*
   Eqn# 6134  yý=a+bx+cxý+dx3+ex4+fx5
 *--------------------------------------------------------------*/
{
  double y;
  y=0.2487113274993664+x*(0.9002141389059909+
    x*(0.9658994322020696+x*(0.007692506399975233+
    x*(-0.0003031672886149768+x*1.121650772708818E-05))));
  return(sqrt(y));
}


/*--------------------------------------------------------------*/
double g31500b(double x)
/*--------------------------------------------------------------*
   Eqn# 4874  y=a+bx+cûxlnx+dx/lnx+e/xý
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=sqrt(x)*log(x);
  x3=x/log(x);
  x4=1.0/(x*x);
  y=3201.553847063860+176.8521250328144*x1
    +436.6821438580977*x2-1883.246271847646*x3
    +1790.467302430289*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31700a(double x)
/*--------------------------------------------------------------*
   Eqn# 4355  y=a+bx+cxý+dxýûx+eexp-x
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x;
  x3=x*x*sqrt(x);
  x4=exp(-x);
  y=0.3578377562139084+1.043214752557405*x1
    -0.001370303280372795*x2-8.509201577454290E-05*x3
    +0.1580591882752711*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31700b(double x)
/*--------------------------------------------------------------*
   Eqn# 4609  y=a+bx+cxýûx+dûxlnx+ex/lnx
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x*sqrt(x);
  x3=sqrt(x)*log(x);
  x4=x/log(x);
  y=2776.473019059725+153.8288087656879*x1
    -0.0001010951827792066*x2+371.9483315265095*x3
    -1620.093039992771*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31800a(double x)
/*--------------------------------------------------------------*
   Eqn# 6132  yý=a+bx+cxý+dx3
 *--------------------------------------------------------------*/
{
  double y;
  y=1.071980333417974+x*(1.360477850187146+
    x*(1.201854443540851+x*-0.008273056193800057));
  return(sqrt(y));
}


/*--------------------------------------------------------------*/
double g31800b(double x)
/*--------------------------------------------------------------*
   Eqn# 4609  y=a+bx+cxýûx+dûxlnx+ex/lnx
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x*sqrt(x);
  x3=sqrt(x)*log(x);
  x4=x/log(x);
  y=5552.947642901565+306.6577058442273*x1
    -0.0002021904025218839*x2+743.8968800012295*x3
    -3240.187020515272*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31900a(double x)
/*--------------------------------------------------------------*
   Eqn# 4382  y=a+bx+cxý+dexpx+eexp-x
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x;
  x3=exp(x);
  x4=exp(-x);
  y=1.948509314122392+1.124895722282613*x1
    -0.005100676601005492*x2+8.905547777852358E-20*x3
    +0.2620286583420982*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g31900b(double x)
/*--------------------------------------------------------------*
   Eqn# 4173  y=a+bx+cxûx+dxý+exýûx
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*sqrt(x);
  x3=x*x;
  x4=x*x*sqrt(x);
  y=28.76729089692713-1.511951568164489*x1
    +0.4217514053565216*x2-0.02633183039679114*x3
    +0.0005859070258838156*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g32000a(double x)
/*--------------------------------------------------------------*
   Eqn# 4341  y=a+bx+cxý+dxýûx+ex3
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x;
  x3=x*x*sqrt(x);
  x4=x*x*x;
  y=3.367449306255549+1.083974300132994*x1
    +0.007668482694768896*x2-0.003614577217487989*x3
    +0.0002675912148501817*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g32000b(double x)
/*--------------------------------------------------------------*
   Eqn# 4755  y=a+bx+cx3+d/lnx+eexp-x
 *--------------------------------------------------------------*/
{
  double y;
  double x1,x2,x3,x4;
  x1=x;
  x2=x*x*x;
  x3=1.0/log(x);
  x4=exp(-x);
  y= -111.6584390050989+1.238144218990064*x1
    -1.739987959934465E-06*x2+379.1717488220770*x3
    -5.512040174539097E+20*x4;
  return(y);
}

/*--------------------------------------------------------------*/
double g4600lo(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/192.8242798917865)/1.748892433468639;
  y=6.966628145396413+65.4192874088499*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4600md(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/1610.269344585912)/2.412647413894015;
  y=11.80584752343706+145.1618675379673*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4600hi(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/2159.088827996704)/2.390534289077195;
  y=14.89281072738765+194.5261398416724*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4700lo(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/216.2009556040568)/1.81202656162805;
  y=6.221403214623654+61.83553855583578*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4700md(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/843.7712566576433)/2.143231970736562;
  y=10.62087345466205+120.3071747968938*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4700hi(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/1308.435220755216)/2.269455130483546;
  y=12.5226863544613+160.3933412357393*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4800lo(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/253.0830911423353)/1.896023727845998;
  y=5.45448266757961+58.6461017619521*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4800md(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/547.1226761430743)/1.946001002988131;
  y=9.179219104642618+105.6311973262788*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4800hi(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/848.3773713172653)/2.154869886473;
  y=10.21004190871538+136.7485496998884*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4900lo(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/206.2626505491258)/1.814962091988218;
  y=3.966946508633288+47.66100216364151*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4900md(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/544.0978143522045)/2.000706807758964;
  y=6.381382418031537+88.54320780795105*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g4900hi(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/1192.45753897097)/2.28873947108425;
  y=9.099751896903934+127.6089430083535*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g41000lo(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/161.725408818154)/1.710574763630868;
  y=2.509991704520077+37.42399134558683*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g41000md(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/525.206855268004)/2.070941812038297;
  y=3.497497088121431+71.24103374128234*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g41000hi(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/2357.682971305255)/2.538559054522592;
  y=7.891852884597797+126.9570676813226*exp(-0.5*n*n);
  return(y);
}
/*--------------------------------------------------------------*/
double g31100lo(double x)
/*--------------------------------------------------------------*
   Eqn# 7203  y=(a+clnx+elnxý)/(1+blnx+dlnxý)
 *--------------------------------------------------------------*/
{
  double y;
  x=log(x);
  y=(1.291826915669442+x*(0.1581477301288331+
    x*0.3560512548608463))/
    (1.0+x*(-0.3816865757039243+x*0.05135364688196362));
  return(y);
}

/*--------------------------------------------------------------*/
double g31100md(double x)
/*--------------------------------------------------------------*
   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
 *--------------------------------------------------------------*/
{
  double y;
  double n;
  n=log(x/461.9583952258151)/2.149631748294819;
  y=0.5145364592257406+53.63085253740166*exp(-0.5*n*n);
  return(y);
}

/*--------------------------------------------------------------*/
double g31100hi(double x)
/*--------------------------------------------------------------*
   Eqn# 7203  y=(a+clnx+elnxý)/(1+blnx+dlnxý)
 *--------------------------------------------------------------*/
{
  double y;
  x=log(x);
  y=(7.934004974027670+x*(-0.2983586873533423+
    x*0.5901343665313362))/
    (1.0+x*(-0.2113457982134224+x*0.01580693446540503));
  return(y);
}

/*--------------------------------------------------------------*/
double g41159lo(double x)
/*--------------------------------------------------------------*
 Uses Extrapolation to estimate moisture content at 1159.
 *--------------------------------------------------------------*/
{ double diff;
  double adjmc;

  diff = (g41000lo(x)-g31100lo(x));
  if (diff < 0.0) diff = (diff * -1);
  adjmc = (g31100lo(x) - ((59.0/60.0)*diff));
  return(adjmc);
}

/*--------------------------------------------------------------*/
double g41159md(double x)
/*--------------------------------------------------------------*
 Uses Extrapolation to estimate moisture content at 1159.
 *--------------------------------------------------------------*/
{ double diff;
  double adjmc;

  diff = (g41000md(x)-g31100md(x));
  if (diff < 0.0) diff = (diff * -1);
  adjmc = (g31100md(x) - ((59.0/60.0)*diff));
 
  return(adjmc);
}


/*--------------------------------------------------------------*/
double g41159hi(double x)
/*--------------------------------------------------------------*
 Uses Extrapolation to estimate moisture content at 1159.
 *--------------------------------------------------------------*/
{ double diff;
  double adjmc;

  diff = (g41000hi(x)-g31100hi(x));
  if (diff < 0.0) diff = (diff * -1);
  adjmc = (g31100hi(x) - ((59.0/60.0)*diff));
 
  return(adjmc);
}



/*----------------------------------------------------------------------------*/
double bdl_ffmc(double ff_ffmc, int hour, int rh)
/*----------------------------------------------------------------------------*
 This function first checks that the program inputs are within the
 appropriate ranges. If any of the input values are out of
 range an error message is printed to the screen and control is
 returned to the operating system.
 NOTE: Appropriate ranges for inputs are defined as follows:
        hour    ( >= 1    and <= 2459  ) (integer)
        RH      ( >= 0    and <= 100   ) (integer)
    FF_FFMC ( >= 17.5 and <= 100.9 ) (one decimal place)
 Where: FF_FFMC = 17.5  corresponds to the lower limit of Van Wagner's
            original diurnal adjustment graphs;
    FF_FFMC = 100.9 corresponds to the theoretical upper limit of
            the current FF_FFMC scale.

 The function then calculates the moisture content at
 1600 (mc1600) based on the std. FF - FFMC value.
 This moisture content is then used in the appropriate equation
 to predict the moisture content at the desired time (adj_mc)
 between 1200 noon on the current day and 1159 the next day.
 Calculations for the various hours are based on the following:
  1200 to 2000 : Equations for every hour (interpolated for minute resolution).
  2001 to  559 : Interpolation between 2000 and 600 (using high RH equation)
  600  to 1100 : Equations for every hour (interpolated for minute resolution)
  1101 to 1159 : Extrapolation using 1100 and 1000 as end-points.

 The function then calculates the new FFMC from the new calculated
 moisture content and returns this time adjusted FFMC value.
*-----------------------------------------------------------------------------*/
{
/*------------- Declare and Initialize Local Variables to 0 ------------------*/
  double mc1600      =0.0 ;
  double adj_mc      =0.0 ;
  double adj_ffmc    =0.0 ;
  double next_hr     =0.0 ;
  double current_hr  =0.0 ;
  int temptime       =0   ;
  double minutes        =0.   ;
  int hrs_since_2000 =0   ;
  int min_since_2000 =0   ;
  double tran_str    =0.0 ;
  double tran_end    =0.0 ;
  double tran_dist   =0.0 ;
  double tran_crve   =0.0 ;
  float  uemcb       =30.0;       /* Lower Equilibrium Moisture Content Bound */
  float  lemcb       =5.0 ;       /* Upper Equilibrium Moisture Content Bound */

/*------------------- Check validity of input data ---------------------------*/
  if (rh <0)
    rh = 90;

  if (hour < 0 || hour > 4859 || rh <0 || rh > 100 || ff_ffmc < 0 || ff_ffmc >= 101.0)
    return(-1.);

  if (hour > 2459) hour-=2400;

  minutes = (double)(hour % 100);                            /*--- Check Time        ---*/
  if (minutes > 59.)
    return(-1.);

  if (ff_ffmc < 17.5) ff_ffmc = 17.5;              /*--- Check for low     ---*/
                                                   /*--- end of FFMC scale ---*/

/*------------------- Calculate moisture content at 1600 ---------------------*/
  mc1600 = (147.2 * (101 - ff_ffmc) / (59.5 + ff_ffmc));

/*------------------- Select the appropriate set of equations ----------------*/

  if (hour >= 600 && hour < 700)
    { 
      if (rh > 87)
        { current_hr    = g4600hi(mc1600);
          next_hr       = g4700hi(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (rh >= 68 && rh <= 87)
        { if (rh > 77)
            { if (hour <= 630) adj_mc = g4600md(mc1600);
              else             adj_mc = g4700hi(mc1600);
            }
      else if (rh >= 58 && rh <= 77)
            { current_hr    = g4600md(mc1600);
              next_hr       = g4700md(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
      else if (rh < 68)
        { if (rh >= 58 && rh <= 77)
            { if (hour <= 630) adj_mc = g4600lo(mc1600);
              else             adj_mc = g4700md(mc1600);
            }
          else if (rh < 58)
            { current_hr    = g4600lo(mc1600);
          next_hr       = g4700lo(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
    }

  else if (hour >= 700 && hour < 800)
    { 
      if (rh > 77)
        { current_hr    = g4700hi(mc1600);
          next_hr       = g4800hi(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (rh >= 58 && rh <= 77)
        { if (rh > 67)
            { if (hour <= 730) adj_mc = g4700md(mc1600);
              else             adj_mc = g4800hi(mc1600);
            }
          else if (rh >= 48 && rh <= 67)
            { current_hr    = g4700md(mc1600);
              next_hr       = g4800md(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
      else if (rh < 58)
        { if (rh >= 48 && rh <= 67)
        { if (hour <= 730) adj_mc = g4700lo(mc1600);
              else             adj_mc = g4800md(mc1600);
            }
          else if (rh < 48)
            { current_hr    = g4700lo(mc1600);
              next_hr       = g4800lo(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
    }

  else if (hour >= 800 && hour < 900)
    { 
      if (rh > 67)
        { current_hr    = g4800hi(mc1600);
          next_hr       = g4900hi(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (rh >= 48 && rh <= 67)
        { if (rh > 62)
            { if (hour <= 830) adj_mc = g4800md(mc1600);
              else             adj_mc = g4900hi(mc1600);
            }
          else if (rh >= 43 && rh <= 62)
            { current_hr    = g4800md(mc1600);
              next_hr       = g4900md(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
      else if (rh < 48)
        { if (rh >= 43 && rh <= 62)
            { if (hour <= 830) adj_mc = g4800lo(mc1600);
              else             adj_mc = g4900md(mc1600);
            }
          else if (rh < 43)
            { current_hr    = g4800lo(mc1600);
              next_hr       = g4900lo(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
    }
    }

  else if (hour >= 900 && hour < 1000)
    { 
      if (rh > 62)
        { current_hr    = g4900hi(mc1600);
          next_hr       = g41000hi(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (rh >= 43 && rh <= 62)
        { if (rh > 57)
            { if (hour <= 930) adj_mc = g4900md(mc1600);
          else             adj_mc = g41000hi(mc1600);
            }
          else if (rh >= 38 && rh <= 57)
            { current_hr    = g4900md(mc1600);
              next_hr       = g41000md(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
      else if (rh < 43)
        { if (rh >= 38 && rh <= 57)
            { if (hour <= 930) adj_mc = g4900lo(mc1600);
              else             adj_mc = g41000md(mc1600);
            }
      else if (rh < 38)
            { current_hr    = g4900lo(mc1600);
              next_hr       = g41000lo(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
            }
        }
    }  

  else if (hour >= 1000 && hour < 1100)
    { if (mc1600 <= 5.0)                                  /*  1000**() */
        { if      (rh > 57)              adj_mc = g41000hi(mc1600);
          else if (rh >= 38 && rh <= 57) adj_mc = g41000md(mc1600);
          else if (rh < 38)              adj_mc = g41000lo(mc1600);
    }
      else if (mc1600 > 5.0 && mc1600 <= 30.0)
        { if (rh > 57)                                /* 1000hi() && 1100hi() */
            { current_hr = g41000hi(mc1600);
              tran_str   = g41000hi(5.0);
              tran_end   = g31100hi(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
          next_hr    = tran_crve;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (rh >=38 && rh <= 57)
        { if (rh > 54.5)                          /* 1000md() && 1100hi() */
        { current_hr = g41000md(mc1600);
          tran_str   = g41000md(5.0);
          tran_end   = g31100hi(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  next_hr    = tran_crve;
                  adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
                }
              else if (rh >= 35.5 && rh <= 54.5)      /* 1000md() && 1100md() */
                { current_hr = g41000md(mc1600);
                  tran_str   = g41000md(5.0);
                  tran_end   = g31100md(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
          next_hr    = tran_crve;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
        }
      else if (rh < 38)
        { if (rh >= 35.5 && rh <= 54.5)           /* 1000lo() && 1100md() */
        { current_hr = g41000lo(mc1600);
          tran_str   = g41000lo(5.0);
          tran_end   = g31100md(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  next_hr    = tran_crve;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
                }
              else if (rh < 35.5)                     /* 1000lo() && 1100lo() */
                { current_hr = g41000lo(mc1600);
                  tran_str   = g41000lo(5.0);
                  tran_end   = g31100lo(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
          next_hr    = tran_crve;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
                }
            }      
        }
      else if (mc1600 > 30.0)
        { if (rh > 57)                                /* 1000hi() && 1100hi() */
            { current_hr    = g41000hi(mc1600);
              next_hr       = g31100hi(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
          else if (rh >= 38 && rh <= 57)
            { if (rh > 54.5)                          /* 1000md() && 1100hi() */ 
                { if (hour <= 1030) adj_mc = g41000md(mc1600);
              else              adj_mc = g31100hi(mc1600);
            }
          else if (rh >= 35.5 && rh <= 54.5)      /* 1000md() && 1100md() */ 
            { current_hr    = g41000md(mc1600);
          next_hr       = g31100md(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
                }
            }
          else if (rh < 38)
            { if (rh >= 35.5 && rh <= 54.5)           /* 1000()lo && 1000()md */ 
                { if (hour <= 1030) adj_mc = g41000lo(mc1600);
              else              adj_mc = g31100md(mc1600);
            }
          else if (rh < 35.5)                     /* 1000lo() && 1100lo() */
            { current_hr    = g41000lo(mc1600);
              next_hr       = g31100lo(mc1600);
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
            }
        }
    }

  else if (hour >= 1100 && hour < 1200)
    { if (mc1600 <= 5.0)                                   /* 1000**() */
        { if      (rh > 57)             adj_mc = g41000hi(mc1600);
          else if (rh >= 38 && rh <=57) adj_mc = g41000md(mc1600);
          else if (rh < 38)             adj_mc = g41000lo(mc1600);
        }
      else if (mc1600 > 5.0 && mc1600 <= 30.0)
    { if (rh > 57)                                /* 1100hi() && 1159hi() */
        { tran_str   = g41000hi(5.0);
              tran_end   = g31100hi(30.0);
              tran_dist  = ((mc1600-lemcb));
              tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
              current_hr = tran_crve;

              tran_end   = g41159hi(30.0);
              tran_dist  = ((mc1600-lemcb));
              tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
              next_hr    = tran_crve;
              adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
            }
      else if (rh >=38 && rh <= 57)
        { if (rh > 54.5)                          /* 1100md() && 1159hi() */
                { tran_str   = g41000md(5.0);
                  tran_end   = g31100hi(30.0);
                  tran_dist  = ((mc1600-lemcb));
                  tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  current_hr = tran_crve;

                  tran_end   = g41159hi(30.0);
                  tran_dist  = ((mc1600-lemcb));
                  tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  next_hr    = tran_crve;
                  adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
                }
          else if (rh >= 35.5 && rh <= 54.5)      /* 1100md() && 1159md() */
                { tran_str   = g41000md(5.0);
                  tran_end   = g31100md(30.0);
                  tran_dist  = ((mc1600-lemcb));
                  tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  current_hr = tran_crve;

                  tran_end   = g41159md(30.0);
                  tran_dist  = ((mc1600-lemcb));
                  tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  next_hr    = tran_crve;
                  adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
                }
        }
      else if (rh < 38)
        { if (rh >= 35.5 && rh <= 54.5)           /* 1000lo() && 1100md() */
                { tran_str   = g41000lo(5.0);
                  tran_end   = g31100md(30.0);
                  tran_dist  = ((mc1600-lemcb));
                  tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  current_hr = tran_crve;
                  
                  if (rh <= 52) tran_end   = g41159md(30.0);
                  else          tran_end   = g41159hi(30.0);
                  tran_dist  = ((mc1600-lemcb));
                  tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
                  next_hr    = tran_crve;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
                }
          else if (rh < 35.5)                     /* 1000lo() && 1100lo() */
        { tran_str   = g41000lo(5.0);
          tran_end   = g31100lo(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
          current_hr = tran_crve;

                  if (rh <  33) tran_end   = g41159lo(30.0);
                  else          tran_end   = g41159md(30.0);
          tran_dist  = ((mc1600-lemcb));
          tran_crve  = tran_str+interp(tran_str,tran_end,tran_dist,(uemcb-lemcb));
          next_hr    = tran_crve;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
        }
        }
    }
      else if (mc1600 > 30.0)
    { if (rh < 35.5)
        { if      (rh < 33.0)          next_hr = g41159lo(mc1600);
          else if (rh >=33 && rh <=52) next_hr = g41159md(mc1600);
          current_hr    = g31100lo(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
              if (adj_mc > current_hr) adj_mc = current_hr;
        }
      else if (rh >=35.5 && rh <=54.5)
        { if      (rh >52)             next_hr = g41159hi(mc1600);
          else if (rh >=33 && rh <=52) next_hr = g41159md(mc1600);
          current_hr    = g31100md(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
        }
      else if (rh > 54.5)
        { next_hr    = g41159hi(mc1600);
          current_hr = g31100hi(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,59);
        }
        }
    }
  else if (hour >= 1200 && hour < 1300)
    {
      if (mc1600 < 21)
        { next_hr       = g31300a(mc1600);
          current_hr    = g31200a(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=21)
        { if      (mc1600  < 22) next_hr = g31300a(mc1600);
          else if (mc1600 >= 22) next_hr = g31300b(mc1600);
          current_hr    = g31200b(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }

  else if (hour >= 1300 && hour < 1400)
    {
      if (mc1600 < 22)
        { next_hr       = g31400a(mc1600);
          current_hr    = g31300a(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=22)
        { if      (mc1600  < 23) next_hr = g31400a(mc1600);
          else if (mc1600 >= 23) next_hr = g31400b(mc1600);
          current_hr    = g31300b(mc1600);
      adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }
  else if (hour >= 1400 && hour < 1500)
    {
      if (mc1600 < 23)
        { next_hr       = g31500a(mc1600);
          current_hr    = g31400a(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=23)
        { next_hr       = g31500b(mc1600);
          current_hr    = g31400b(mc1600);
      adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }

  else if (hour >= 1500 && hour < 1600)
    {
      if (mc1600 < 23)
        { next_hr       = mc1600;
          current_hr    = g31500a(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=23)
        { next_hr       = mc1600;
      current_hr    = g31500b(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }

  else if (hour >= 1600 && hour < 1700)
    {
      if (mc1600 < 40)
        { next_hr       = g31700a(mc1600);
          current_hr    = mc1600;
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
         }
      else if (mc1600 >=40)
     { next_hr       = g31700b(mc1600);
           current_hr    = mc1600;
           adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
         }
    }
  else if (hour >= 1700 && hour < 1800)
    {
      if (mc1600 < 40)
        { next_hr       = g31800a(mc1600);
          current_hr    = g31700a(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=40)
    { next_hr       = g31800b(mc1600);
          current_hr    = g31700b(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }

  else if (hour >= 1800 && hour < 1900)
    {
      if (mc1600 < 40)
        { next_hr       = g31900a(mc1600);
          current_hr    = g31800a(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=40)
        { if      (mc1600  < 42) next_hr = g31900a(mc1600);
          else if (mc1600 >= 42) next_hr = g31900b(mc1600);
          current_hr    = g31800b(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }

  else if (hour >= 1900 && hour < 2000)
    {
      if (mc1600 < 42)
        { next_hr       = g32000a(mc1600);
          current_hr    = g31900a(mc1600);
      adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
      else if (mc1600 >=42)
        { if      (mc1600  < 49) next_hr = g32000a(mc1600);
          else if (mc1600 >= 49) next_hr = g32000b(mc1600);
          current_hr    = g31900b(mc1600);
          adj_mc = current_hr + interp(current_hr,next_hr,minutes,60.);
        }
    }

  else if (hour == 2000)
    {
      if      (mc1600 < 49) adj_mc = g32000a(mc1600);
      else if (mc1600 >=49) adj_mc = g32000b(mc1600);
    }

  else if (hour > 2000 || hour < 600)        /*--------- Extrapolate ---------*/
    {
      if      (mc1600 < 49) current_hr = g32000a(mc1600);
      else if (mc1600 >=49) current_hr = g32000b(mc1600);

      next_hr = g4600hi(mc1600);      /* Uses High RH Curve as Default */

      if (hour < 2000)
        { temptime       = hour + 2400;
          hrs_since_2000 = ((temptime - 2000) - (int)minutes) /100;
      min_since_2000 = (hrs_since_2000 * 60) + (int)minutes;
      adj_mc  = current_hr + interp(current_hr,next_hr,min_since_2000,600.);
        }
      else
    { hrs_since_2000 = ((hour - 2000) - (int)minutes) /100;
      min_since_2000 = (hrs_since_2000 * 60) + (int)minutes;
      adj_mc  = current_hr + interp(current_hr,next_hr,min_since_2000,600.);
    }
    }
  else
    { return(-1);
    }

  adj_ffmc =  (59.5 * (250 - adj_mc) / (147.2 + adj_mc));
  if (adj_ffmc < 0)     adj_ffmc = 0;
  if (adj_ffmc > 101.0) adj_ffmc = 101.0;
 
  return(adj_ffmc);
}
