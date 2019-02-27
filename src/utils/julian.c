#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* This routine converts the time in decimal hours to military time HHMM (e.g. 1200)*/
   
int mtime(double dhour)
{
double x, n;
int hour;

    while (dhour>=24.) dhour -= 24.;
    while (dhour<0.) dhour += 24.;
    
    x = modf(dhour,&n); 
    hour = (int)(100*n+60*x);
//    printf("... %d %lf\n", hour, dhour);

    return(hour);
}
/* Two definitive routines to convert to and from julian day.
   Note that months start with Jan=1 */

int Date2Julian(int month, int day, int year)
{
struct tm *ptTime;
time_t lTime;
int julian;
int days[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

   if (year == 0)
   {
      time(&lTime);
      ptTime = localtime(&lTime);
      year = 1900+ptTime->tm_year;
   }

   julian = days[month-1]+day;

   if ((year%4) == 0 && month > 2)
     julian++;

   return(julian);
}
/***************************************************/
void Julian2Date(int julian, int *year, int *month, int *day)
{
struct tm *ptTime;
time_t lTime;
char *months[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
int i;
int days[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

   *month=0; *day=0;
   if (*year == 0)
   {
      time(&lTime);
      ptTime = localtime(&lTime);
      *year = 1900+ptTime->tm_year;
   }

   if (*year%4 == 0)
      for (i=2; i<12; i++)
         days[i]++;

   if (julian > days[11])
      *month = 12;
   else
      while (julian > days[*month] && *month < 12) (*month)++;

   *day = julian - days[*month-1];
}
