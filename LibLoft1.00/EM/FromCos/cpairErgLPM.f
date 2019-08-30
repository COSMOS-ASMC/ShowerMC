!     test cpairErgLPM
!
!      real*8 ee, rho, eg
!      integer i
!      read(*,*) eg, rho
!      eg = eg/1.e9      ! GeV
!      rho = rho * 1.e3   ! kg/m3
!      do i = 1, 10000
!         call cpairErgLPM(eg, rho, ee)
!         write(*, *) sngl(ee/eg)
!      enddo
!      end
!      
!       sampling  Pair electron energy
!       when the LPM effect is considered.
!       this may be used eg > 10^17 eV 
!

      subroutine cpairErgLPM(eg, rho, ee)
      implicit none
#include "Zmass.h"
      real*8 eg  !input.  gamma energy in GeV
      real*8 rho !input.  air density in kg/m3
      real*8 ee  !output. sampled electron or positron energy in GeV
!
      real*8  u, v, temp, s, w, u1, v1, w1, x
      real*8 cmigdG, cmigdPsi, cmigdSforPair, cmigdGzai

      real*8 a1, a2,  suma, a1prob
      parameter ( a1 = 2./3., a2 = 1./9.,
     *           suma = a1 + a2, a1prob = a1/suma)
!
      do while (.true.)
         call rndc(u)
         if(u .lt. a1prob) then
!             sampling by  dv
            call rndc(v)
!            rejection
!                get s
            s = cmigdSforPair(eg, rho, v)
            temp = cmigdGzai(s) * (cmigdG(s) +cmigdPsi(s))/2.
            call rndc(u)
            if(u .lt. temp) then
!                  accepted
               ee =max(min(v, 1.d0-masele/eg),masele/eg) * eg
               goto 100
            endif
         else
!              sampling by 12(v-0.5)dv
            call rndc(u)
            call rndc(v)
            call rndc(w)
!             get one which is most far from 1/2
            u1 = abs(u - 0.5)
            v1 = abs(v - 0.5)
            w1 = abs(w - 0.5)
            x = max(u1, v1, w1)
            if(x .eq. u1) then
               v = u
            elseif(x .eq. w1) then
               v = w
            endif
!            rejection, get s
            s = cmigdSforPair(eg, rho, v)
            temp = cmigdGzai(s) * cmigdPsi(s)
            call rndc(u)
            if(u .lt. temp) then
!               accepted
               ee =max(min(v, 1.d0-masele/eg), masele/eg) * eg
               goto 100
            endif
         endif
      enddo
 100  continue
      end
