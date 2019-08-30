      subroutine cmuPrsmpP(Emu, prob, path)
      implicit none
#include "Zcmuint.h"

      real*8 Emu             ! input.  muon total energy in GeV
      real*8 prob            ! output. muon pair creation prob. /X0
      real*8 path            ! output. sampled path in r%l
                             ! if Emu < media.cnst.muPrEmin,  prob=0 
                             !  and path becomes big

      real*8 u, ale

      if(Emu .le. muPrEmin) then
         prob = 0.
      elseif(Emu .le. muPrEmax1) then
         ale = log10(Emu)
         call kintp3(MuPrTX, 1, muPrTXT, muPrLEmin,
     *   muPrdETX, ale, prob)
      else
!          prob is const for  v > vmin; Emu > muPrEmax1         
         prob = MuPrTX(muPrTXT)
      endif
      if(prob .gt. 0.) then
         call rndc(u)
         path =- log(u)/prob
      else
         path = 1.d30
      endif
      end
      
      subroutine cmuPrsmpE(Emu, Epair)
      implicit none
#include "Zcmuint.h"
      real*8 Emu           ! input. muon total energy in GeV
      real*8 Epair         ! output. sampled energy loss of muon
                           !        (Epair/Emu > media.cnst.muPrVmin)

      real*8  a, b, v

      a = mupa*Emu**(-0.40) + mura
      b = mupb*Emu**muqb  +  0.92

      call cmuPrsmp0(a, b, muPrVmin, v)
      Epair = Emu*v
      end


!     **********************************
      subroutine cmuPrsmp0(a, b, vc, v)
!     **********************************
      implicit none
!       sample v from  dv/v**(2-b)/(1+av**b)**2
!       This is a good approximation for pair spectrum
!       from muon 
      real*8 a     ! input.
      real*8 b     ! input.
      real*8 vc    ! input.  minimum of v
      real*8 v     ! output. sampled v


      real*8 u, temp


      logical ok

      ok = .false.
      do while (.not. ok)
         if(b .ne. 1.0d0) then
            temp = vc**(b-1)
            call rndc(u)
            v = ((1-temp)*u + temp)**(1./(b-1))
            call rndc(u)
            ok =(u .lt. ((1.+a*temp*vc)/(1.+a*v**b))**2) 
         else
            call rndc(u)
            v = vc**u
            call rndc(u)
            ok =(u .lt. ((1.+a*vc)/(1.+a*v))**2) 
         endif
      enddo

      end
