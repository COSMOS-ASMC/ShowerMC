      subroutine epXrayp(media, Exin, m1, m2, p, path)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"


       type(epmedia):: media     ! input.  media
      real*8 Exin               ! input. X/Gamma energy in GeV  (mostly keV to MeV). <100.
      integer m1, m2            ! input. 1~5.    1 --> coh. scatt.
                                !                2 --> + incoh. scatt.
                                !                3 --> + p.e  
                                !                4 --> + pair cre. by nuc.
                                !                5 --> + pari cre. by elec.
                                !             are specified
      real*8 p(m2)              ! output. probability ( number of
                                ! occurence ) of m-th process per r.l
      real*8 path(m2)           ! output. sampled path in r.l for m-th process.

                              
      real*8 u
      integer icon, i
      real Ex, xsec(7)
      Ex=Exin
      if(media%xcom%size .gt. 0) then
         call cGetXXsec(Ex, media%xcom%tab, media%xcom%size, 
     *    m1, m2,  xsec, icon)
!          xsec(m) is 1/(g/cm2).
         do i = m1, m2
            p(i) = xsec(i)*media%X0g
            if( p(i) .gt.  0.) then
               call rndc(u)
               path(i) = -log(u)/p(i)
            else
               path(i) = 1.e35
            endif
         enddo      
      else
         write(0,*) 'size=0 in epXrayp%f for media=', media%name
      endif
      end
