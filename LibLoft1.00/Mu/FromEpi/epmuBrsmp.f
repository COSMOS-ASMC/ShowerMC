      subroutine epmuBrsmpP(media, Emu, prob, path)
      implicit none
#include "Zmedia.h"

       type(epmedia):: media    !input. media
      real*8 Emu             ! input.  muon total energy in GeV
      real*8 prob            ! output. muon brems creation prob. /X0
      real*8 path            ! output. sampled path in r.l
                             ! if Emu < media.cnst.muBrEmin,  prob=0 
                             ! and path becomes big

      real*8 u, ale

      if(Emu .le. media%cnst%muBrEmin) then
         prob = 0.
      elseif(Emu .le. media%cnst%muBrEmax1) then
         ale = log10(Emu)
         call kintp3(media%tbl%MuBrTX, 
     *   1,  media%cnst%muBrTXT, media%cnst%muBrLEmin,
     *   media%cnst%muBrdETX, ale, prob)
      else
!          prob is const for  v > vmin; Emu > muBrEmax1         
         prob = media%tbl%MuBrTX(media%cnst%muBrTXT)
      endif
      if(prob .gt. 0.) then
         call rndc(u)
         path =- log(u)/prob
      else
         path = 1.d30
      endif
      end
      
      subroutine epmuBrsmpE(media, Emu, Eg)
      implicit none
#include "Zmass.h"
#include "Zmedia.h"
       type(epmedia):: media  ! input.  media
      real*8 Emu           ! input. muon total energy in GeV
      real*8 Eg            ! output. sampled energy loss of muon
!
!         brems v=Eg/Emu distribution:
!          dv epmuBremLogf*(4(1 - v)/3 + v)/v
!
      real*8 vc, vmx, epmuvmax2, term1, term2, u1, u2, x
      real*8 delta, logf, func, u, v

      vc = media%cnst%muBrVmin
      vmx = epmuvmax2(media, Emu)
      term1  =4.d0/3.d0 * (log(vmx/vc) - (vmx-vc))
      term2  =  (vmx-vc)*(vmx+vc)/2.d0
!          loop for rejection
 100  continue
      call rndc(u)
      if(u .le.  term1/(term1+term2)) then
!         (1/x -1)dx
         do while (.true.)
!              average number of trials is 1.0xx; xx depends on vc
            call rndc(u)
            x = vc**u
            call rndc(u)
            if( u .lt. (1.0-x)) then
               if(x .lt. vmx) goto 10
            endif
         enddo
 10      continue
      else
!          x dx in [vc:vmx]
         do while (.true.)
            call rndc(u1)
            call rndc(u2)
            x = max(u1,u2)
            if(x .gt. vc .and. x .lt. vmx) goto 20
         enddo
 20      continue
      endif
      v = x
!         final rejection by epmuBremLogf.
!         we use specially made epmuBremLogf here

!         mim. momentum transfer in unit of Mmu
      delta = (masmu/Emu) * v /(1.d0-v)/2
      logf =  media%mu%Akm / (1.d0 + media%mu%Akm2*delta)
!      if(media.Zeff .gt.  10.) then
      if(media%Z .gt.  10.) then
!         logf = logf * 2./3.d0 /media.Zeff3
         logf = logf * 2./3.d0 /media%Z**0.3333
      endif
      func = log(logf)  ! = epmuBremLogf
      call rndc(u)
      if(u .gt. func/media%mu%logf0)  goto 100
      Eg = v * Emu
      end

      




