!
!             dE/dx fluctuation for heavy ions
!
      subroutine epdedxflhv(media, aPtcl, sigmacol2, sigmachg2)
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"

       type(epmedia):: media       ! input. media

       type(ptcl):: aPtcl        ! input. a particle
      real*8 sigmacol2          ! output.  square of Gaussian simga
                                !   (GeV^2/(g/cm^2)).
                                !   VERY IMPORTANT. You have to
                                !  multiply the length (in g/cm^2) to 
                                !  sigma2 and take its root to get the
                                !  sigma.  
      real*8 sigmachg2          ! output.  same due to charge attachment

      real*8  gamma,  temp
      real*8  epdedxqeff
      real*8  qeff, coef

      character*80 msg
      data coef/0.307d-3/   !  GeV /(g/cm^2)

      if(aPtcl%code .ne. kgnuc) then
         write(msg, *) ' ptcl code=',aPtcl%code,
     *     ' is not a heavy particle in epdedxflhv1'
         call cerrorMsg(msg, 0)
      endif
      gamma = aPtcl%fm%p(4)/aPtcl%mass
      qeff = epdedxqeff(aPtcl)
!cc      temp =  coef*qeff**2 * masele* media.ZbyAeff 
!             qeff gives too small fluctuation
      temp =  coef* aPtcl%charge**2 * masele* media%ZbyAeff 
      sigmacol2 = temp * (1.+ gamma**2)/2.
      sigmachg2 = temp * 1.25*(1.- qeff/aPtcl%charge)
      end
