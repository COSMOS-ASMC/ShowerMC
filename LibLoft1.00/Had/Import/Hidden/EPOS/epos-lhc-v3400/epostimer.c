#include <sys/times.h>
#include <sys/timeb.h>
#include <stdio.h>

void timer_(long *etime){
struct timeb tr;
struct tms tu;
ftime (&tr);
times (&tu);
etime[0]=tu.tms_utime;
etime[1]=tu.tms_stime;
etime[2]=tr.time;
etime[3]=(long) tr.millitm;
etime[4]=(long) tr.timezone;
}
