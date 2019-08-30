      subroutine cmuNsmpP(Emu, prob, path)
      implicit none
#include "Zcmuint.h"

      real*8 Emu             ! input.  muon total energy in GeV
      real*8 prob            ! output. muon nuc. int prob. /X0
      real*8 path            ! output. sampled path in r%l
                             ! if Emu < media.cnst.muNEmin,  prob=0 
                             ! and path becomes big

      real*8 u, ale

      if(Emu .le. muNEmin) then
         prob = 0.
      elseif(Emu .le. muNEmax1) then
         ale = log10(Emu)
         call kintp3(MuNTX, 1, muNTXT, muNLEmin,
     *   muNdETX, ale, prob)
      else
!        small power dependence
         prob = MuNTX(muNTXT)*
     *      (Emu/muNEmax1)**muNpwtx
      endif
      if(prob .gt. 0.) then
         call rndc(u)
         path =- log(u)/prob
      else
         path = 1.d30
      endif
      end
      
      subroutine cmuNsmpE(Emu, Et)
      implicit none
#include "Zcmuint.h"
      real*8 Emu           ! input. muon total energy in GeV
      real*8 Et            ! output. sampled energy transfer
!

      real*8  ale, u, uu,   v


      real*8 a
      real*8 error

      data a/0.02d0/   ! to change this, you must change creation part



      ale = log10(min(Emu,  muNEmax))


      call rndc(u)
      uu = a*u/(1+a-u)  ! uniform in this variable
         call kpolintp2(0.d0, 1, muNdU,
     *        muNLEmin, 1, muNdE,
     *        MuNTbl,  muNUsize,  muNUsize, muNEsize,
     *        5, 3,   uu, ale, v, error)
      Et = 10.d0**v * Emu
      end

