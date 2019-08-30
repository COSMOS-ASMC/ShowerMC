!     test cbremErgLPM
!
!      real*8 ee, rho, eg
!      integer i
!      read(*,*) ee, rho
!      ee = ee/1.e9      ! GeV
!      rho = rho * 1.e3   ! kg/m3
!      do i = 1, 10000
!         call cbremErgLPM(ee, rho, eg)
!         write(*, *) sngl(eg/ee)
!      enddo
!      end
!      
!       sampling  Bremsstrahlung gamma ray energy.
!       when the LPM effect is considered.
!       this may be used for ee > 10^15 eV.
!
!     In unified version (Cosmos9, Epics10.)  constants inside
!     is calculated for each medium, and rho takes into account
!     air density change  or the rhoc parameter,
!     next  can be used other than  Air.
      subroutine cbremErgLPM(ee, rho, eg)
      implicit none
#include "Zmass.h"
      real*8 ee  !input.  electron energy in GeV
      real*8 rho !input.  air density in kg/m3
      real*8 eg  !output.  sampled gamma ray energy. GeV
!
      real*8  u, v, temp, s
      real*8 cmigdG, cmigdPsi, cmigdSforBrem, cmigdGzai

      common /cmigdc/ s1, logs1, vm, logvm, smpa1, smpa2, suma,
     *   a2prob, X0
      real*8 vm, logvm, smpa1, smpa2, suma, a2prob, X0
      real*8 s1  ! /1.1185e-4/  ! (Z^0.333/183)^2 ; Z=7.25
      real*8 logs1   !  /-9.0983/  ! ln(s1)

!      real*8 a1, a2, vm, logvm, suma, a2prob
!      parameter (vm = 1.d-4, logvm = -4*2.3025851, 
!     *           a1 = 0.5, a2 = -4.d0*logvm/3.,
!     *           suma = a1 + a2, a2prob = a2/suma)
!
      do while (.true.)
         call rndc(u)
         if(u .lt. a2prob) then
!             sampling by  1/v
            call rndc(u)
            v = exp(u * logvm)  ! fractional energy
!            rejection
            call rndc(u)
            if(u .lt. 1.-v) then
!                further check; get s
               s = cmigdSforBrem(ee, rho, v)
               if(u .lt. (1.-v) * cmigdGzai(s)*cmigdPsi(s)) then
!                  accepted
                  eg =min( v, 1.d0- masele/ee) * ee
                  goto 100
               endif
            endif
         else
!              sampling by 2vdv
            call rndc(u)
            call rndc(v)
            v = max(u, v)
!            rejection, get s
            s = cmigdSforBrem(ee, rho, v)
            temp = cmigdGzai(s) * (cmigdG(s) + 2* cmigdPsi(s))/3.
            call rndc(u)
            if(u .lt. temp) then
!               accepted
               eg = min(v, 1.d0-masele/ee) * ee
               goto 100
            endif
         endif
      enddo
 100  continue
      end
