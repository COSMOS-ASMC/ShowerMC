      subroutine cmuBrsmpP(Emu, prob, path)
      implicit none
#include "Zcmuint.h"

      real*8 Emu             ! input.  muon total energy in GeV
      real*8 prob            ! output. muon brems creation prob. /X0
      real*8 path            ! output. sampled path in r%l
                             ! if Emu < muBrEmin,  prob=0 
                             ! and path becomes big

      real*8 u, ale

      if(Emu .le. muBrEmin) then
         prob = 0.
      elseif(Emu .le. muBrEmax1) then
         ale = log10(Emu)
         call kintp3(MuBrTX, 1,  muBrTXT, muBrLEmin,
     *    muBrdETX, ale, prob)
      else
!          prob is const for  v > vmin; Emu > muBrEmax1         
         prob = MuBrTX(muBrTXT)
      endif
      if(prob .gt. 0.) then
         call rndc(u)
         path =- log(u)/prob
      else
         path = 1.d30
      endif
      end
      
      subroutine cmuBrsmpE(Emu, Eg)
      implicit none
#include "Zcmuint.h"
#include "Zmass.h"

      real*8 Emu           ! input. muon total energy in GeV
      real*8 Eg            ! output. sampled energy loss of muon
!
!         brems v=Eg/Emu distribution:
!          dv cmuBremLogf*(4(1 - v)/3 + v)/v
!
      real*8 vc, vmx, cmuvmax2, term1, term2, u1, u2, x
      real*8 delta, logf, func, u, v

      vc = muBrVmin
      vmx = cmuvmax2(Emu)
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
!         final rejection by cmuBremLogf.
!         we use specially made cmuBremLogf here

!         mim. momentum transfer in unit of Mmu
      delta = (masmu/Emu) * v /(1.d0-v)/2
      logf =  muAkm / (1.d0 + muAkm2*delta)
      if(Zeff .gt.  10.) then
         logf = logf * 2./3.d0 /Zeff3
      endif
      func = log(logf)  ! = cmuBremLogf
      call rndc(u)
      if(u .gt. func/mulogf0)  goto 100
      Eg = v * Emu
      end
      




