!           ***************************************************
!              cgetLoginN:  get login name ( <= 8 characters)
!           *************************************************** 
!        +U77 compliler option is needed
!        for the HP compiler.
!        
#undef GENERIC

#if defined DECALPHA
#define GENERIC
#elif defined PCLinux
#define GENERIC
#elif defined PCLinuxIFC64
#define GENERIC
#elif defined (PCLinuxIFC)  || defined (KEKA) || defined (KEKB) 
#define GENERIC
#elif defined MACOSX
#define GENERIC
#elif defined CF_AlphaLinux
#define GENERIC
#endif


#if defined NEXT486
!           ******************** next Absoft  fortran ***************
!           the same one as sun/hp may be usable instead of this.
!
       subroutine cgetLoginN(userid)
        character*(*) userid
 	integer getlogin,  intlogin
        pointer(pu, intlogin)
 	pu=getlogin()
         call copyCstring(intlogin, userid)
       end
       subroutine copyCstring(Cstring, target)
         implicit none
         integer i
         character Cstring(*), target*(*)
         
         target=' '
         do i=1, len(target)
         	if(Cstring(i) == char(0) ) exit
         	target(i:i)=Cstring(i)
         enddo
       end
#elif defined GENERIC
      subroutine cgetLoginN(userid)
      implicit none
      character*(*) userid
      character*16 user
      character*1 NULL
      integer leng, kgetenv
      NULL = char(0)
      leng = kgetenv("USER"//NULL, user)
      userid = user(1:leng)
      end
#else
!     **************************** sun4/hp ***************************
!       !!!!  for HP you need +U77 compiler option
!
       subroutine cgetLoginN(userid)
!
       character*(*)  userid
       call getlog(userid)
       end
#endif
#ifdef IBMAIX
        subroutine getlog(user)
        character user*(*)
        call getenv('USER',user)
        end
#endif
!
!       easier interface to kgetenv
!
      integer function kgetenv2(envname, envresult)

      character*(*)  envname   ! input.  environment variable name
      character*(*)  envresult ! output. value of hte env. variable.
 
      integer leng, kgetenv, klena
      character*1  NULL
      character*128 path
      NULL = char(0)
      leng = klena(envname)
      if( leng .gt. 0 ) then
         leng = kgetenv( envname(1:leng)//NULL, envresult)
      else
         leng = 0
         envresult = ' '
      endif
      envresult=envresult(1:leng)//" "
      kgetenv2 = leng
      end
      subroutine cqLibVersion(libv)
      implicit none
#include "Zmanagerp.h"
      character*8::libv ! output  cosmos version such 7.58; left justified
      character*64::LIBLOFT
      character*128:: filen
      integer kgetenv2, icon

      if( kgetenv2("LIBLOFT", LIBLOFT) == 0 ) then
         write(0,*) "Environmental variable "
         write(0,*) "LIBLOFT cannot be obtained in cqLibVersion"
         stop
      endif
      filen = trim(LIBLOFT)//"/Version/version" 
      call copenf(TempDev, filen, icon)
      read(TempDev, '(a)') libv
      close(TempDev)
      end
