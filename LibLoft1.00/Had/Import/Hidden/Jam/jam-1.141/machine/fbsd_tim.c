/* $Id */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/errno.h>
#include <sys/times.h>
#include <string.h>
/* #include <sys/types.h> */


void dtime_(ar)
double ar[5];
{
struct tms timar;
clock_t ct;

  ct=times(&timar);
  ar[0]=(double)ct/CLK_TCK;          /* total time */
  ar[1]=timar.tms_utime/CLK_TCK;     /* user CPU time */
  ar[2]=timar.tms_stime/CLK_TCK;     /* system CPU time */
  ar[3]=timar.tms_cutime/CLK_TCK;    /* system CPU time of child processes */
  ar[4]=timar.tms_cstime/CLK_TCK;    /* system CPU time of child processes */

}

int idtime_()
{
	clock_t	clk;
	clk = clock();
/*	printf("CPU %d sec\n",clk/CLOCKS_PER_SEC) */
	return(clk/CLOCKS_PER_SEC);
}



 
void dateb_(today,ar)

char *today;
int ar[7];
{
time_t tt;
struct tm *tm_buf;

if ((tt = time(NULL)) == -1) {
     fprintf(stderr, "errono:%d %s\n",errno,
       (errno > 0 && errno < sys_nerr ? sys_errlist[errno] : ""));
      exit(1);
}

  tm_buf = localtime(&tt);
  strcpy(today,ctime(&tt));

         ar[0]=tm_buf->tm_wday;
         ar[1]=tm_buf->tm_hour;
         ar[2]=tm_buf->tm_min;
         ar[3]=tm_buf->tm_sec;
         ar[4]=tm_buf->tm_year;
         ar[5]=tm_buf->tm_mon + 1;
         ar[6]=tm_buf->tm_mday;

/*
   printf("%s\n",today);

  fprintf(stderr,"%d:%d:%d  %d-%d-%d\n",
         ar[1],ar[2],ar[3],ar[4],ar[5],ar[6]);

  fprintf(stderr,"%d:%d:%d  %d-%d-%d\n",
         tm_buf->tm_hour,
         tm_buf->tm_min,
         tm_buf->tm_sec,
         1900 + tm_buf->tm_year,
         tm_buf->tm_mon + 1,
         tm_buf->tm_mday);

  fprintf(stderr,"%d/%d/%d %d:%d:%d\n",
         1900 + tm_buf->tm_year,
         tm_buf->tm_mon + 1,
         tm_buf->tm_mday,
         tm_buf->tm_hour,
         tm_buf->tm_min,
         tm_buf->tm_sec);
*/

}


