#include "stdio.h"
#include "stdlib.h"
#include <unistd.h>
#include "string.h"

#if defined NEXT486
int kgetenv(char *env,  char *name){
#elif defined MACOSX
int KGETENV(char *env,  char *name){
#elif defined PCLinux
int KGETENV(char *env,  char *name){ 
#elif defined (IBMAIX) ||  defined (KEKA) || defined (KEKB)
int kgetenv(char *env,  char *name){
#else
int kgetenv_(char *env,  char *name){
#endif
/*    env: input.  Environmental string to be examined.
*                 must have NULL at the end.
*     name: output. Environmental variable's value is put
*                   must have enough length (NULL is not put)
*     function value: length of the name.
*     If env is not defnied.  0 is returned.
*  The fortran routine must reset the string like
*     integer leng
*     character*1  NULL 
*     character*128 path 
*     NULL = char(0)
*     leng = kgetenv("COSMOSTOP"//NULL, path)
*     path = path(1:leng)
*/
  char *envloc;
/*  printf("%s is input\n", env); */
  if((envloc=getenv(env)) == NULL )  {
/*    printf("not found\n"); */
    return 0;
  }
  else {
/*    printf("%s\n", envloc); */
    strcpy(name, envloc);
    return  strlen(envloc);
  }
}
/*   
 *     get unix process number of this program run 
 *     Usage:  
 *     integer procnumber, kgetpid, somedummy
 *     procnumber = kgetpid(somedummy)
*/
#if defined (NEXT486) || defined (IBMAIX) || defined (KEKA) || defined (KEKB)
 int kgetpid(int dummy)
#elif defined (MACOSX) || defined (PCLinux)
   int KGETPID(int dummy)
#else
   int kgetpid_(int dummy)
#endif
   {
     return( (int)getpid());
   }
/*  
 *       get standard unix time
 *      Usage:
 *      integer  inttime, somedummy
 *      inttime = kgettime(domedummy) 
 */
#include <time.h>
#if defined (NEXT486) || defined (IBMAIX) || defined (KEKA)  || defined (KEKB)
 int kgettime(int dummy)
#elif defined (MACOSX) || defined (PCLinux)
   int KGETTIME(int dummy)
#else
   int kgettime_(int dummy)
#endif
   {
     time_t  longtime;     
     longtime = time(NULL);
     return  (int) longtime;
   }

/* 
 *   get present YYMMDDHHMMSS 
 *  Usage:
 *     integer dummy
 *     character*12  yymmdd
 *     dummy = kgetnow(yymmdd)
 *     yymmdd dose not include "\n"
 *     return value is 12 which is the length of yymmdd
 */

#if defined (NEXT486) || defined (IBMAIX) || defined (KEKA)  || defined (KEKB)
 int kgetnow(char *now)
#elif defined (MACOSX) || defined (PCLinux)
   int KGETNOW(char *now)
#else
   int kgetnow_(char *now)
#endif
   {
     struct tm *stime; 
     time_t  longtime;     
     char nowtemp[13];
     longtime = time(NULL);
     stime = localtime( &longtime );
     sprintf( nowtemp,
       "%-.2d%-.2d%-.2d%-.2d%-.2d%-.2d", (stime->tm_year-100),
       stime->tm_mon+1, stime->tm_mday,
       stime->tm_hour, stime->tm_min,stime->tm_sec);
     strncpy(now, nowtemp, 12);
     return 12;
   }
