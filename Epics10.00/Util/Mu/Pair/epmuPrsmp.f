      subroutine epmuPrsmpP(media, Emu, prob, path)
      implicit none
#include "Zmedia.h"

       type(epmedia):: media    !input. media
      real*8 Emu             ! input.  muon total energy in GeV
      real*8 prob            ! output. muon pair creation prob. /X0
      real*8 path            ! output. sampled path in r.l
                             ! if Emu < media.cnst.muPrEmin,  prob=0 
                             !  and path becomes big

      real*8 u, ale

      if(Emu .le. media%cnst%muPrEmin) then
         prob = 0.
      elseif(Emu .le. media%cnst%muPrEmax1) then
         ale = log10(Emu)
         call kintp3(media%tbl%MuPrTX, 
     *   1,  media%cnst%muPrTXT, media%cnst%muPrLEmin,
     *   media%cnst%muPrdETX, ale, prob)
      else
!          prob is const for  v > vmin; Emu > muPrEmax1         
         prob = media%tbl%MuPrTX(media%cnst%muPrTXT)
      endif
      if(prob .gt. 0.) then
         call rndc(u)
         path =- log(u)/prob
      else
         path = 1.d30
      endif
      end
      
      subroutine epmuPrsmpE(media, Emu, Epair)
      implicit none
#include "Zmedia.h"
       type(epmedia):: media  ! input.  media
      real*8 Emu           ! input. muon total energy in GeV
      real*8 Epair         ! output. sampled energy loss of muon
                           !        (Epair/Emu > media.cnst.muPrVmin)


      real*8  a, b, v

      a = media%mu%pa*Emu**(-0.40) + media%mu%ra
      b = media%mu%pb*Emu**media%mu%qb  +  0.92

      call epmuPrsmp0(a, b, media%cnst%muPrVmin, v)
      Epair = Emu*v
      end


!     **********************************
      subroutine epmuPrsmp0(a, b, vc, v)
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


