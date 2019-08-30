      subroutine epBrCSampP(media,  Ee, prob)
      use modEMcontrol
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmedia.h"
#include "Zmass.h"
! #include "ZepTrackp.h"
      type(epmedia)::  media  ! input
      real*8 Ee  ! input. Electron energy in GeV
      real*8 prob  !  output. prob. of Brems /  r.l
      
      real*8 vc, Ek
      Ek = Ee - masele
!      vc = media.cnst.BremEgmin/(Ee-masele)
      vc = media%cnst%BremEgmin
      prob = 4. * ar02 * ( ( log(1.d0/vc) -
     *   (1.d0-vc)) * media%cScrMain +
     +   (1.d0-vc)* (1.d0+vc)/2 * media%cScrC1 )
     *   * media%mbtoPX0

      if( HowNormBrems == -1 ) then
         ! nothing to do
      else
         if( HowNormBrems == 1 ) then
            prob = prob/media%cnst%NormS
         elseif( HowNormBrems == 0 ) then
          !   prob = prob/media.cnst.NormCS   no need normcs=1
         else
            write(0,*) ' HowNormBrems= ', HowNormBrems
            write(0,*) ' invalid in epBrCSampP '
            stop
         endif
      endif
      end

!
!     ************
      subroutine epBrCSampE(media, Ee, Eg)
!     ************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media  ! input
      real*8 Ee  ! input. Electron energy in GeV
      real*8 Eg  !  output. sampled Eg

      real*8 vc, x, u, u1, u2, term1, term2

!      vc = media.cnst.BremEgmin/(Ee-masele)
      vc = media%cnst%BremEgmin


      term1  = media%cScrMain * (log(1.d0/vc) - (1.d0-vc))
      term2  = media%cScrC1 * (1.0-vc)*(1.0+vc)/2.d0
      call rndc(u)
      if(u .le.  term1/(term1+term2)) then
!         (1/x -1)dx
         do while (.true.)
!              average number of trials is 1.0xx; xx depends on vc
            call rndc(u)
            x = vc**u
            call rndc(u)
            if( u .lt. (1.0-x)) goto 10
         enddo
 10      continue
      else
!          x dx in [vc:1]
         do while (.true.)
            call rndc(u1)
            call rndc(u2)
            x = max(u1,u2)
            if(x .gt. vc ) goto 20
         enddo
 20      continue
      endif
      Eg = (Ee-masele) * x
      end
