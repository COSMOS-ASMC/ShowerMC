!#include "cblkdedx.h"
!      implicit none
!      real*8 e, dedt, ek
!      integer i
!      call cdedxEleci(225.d-6,  .true.)
!      do i=1, 300
!        ek=1.e-9*10.**( (i-1) /20.) 
!        e = ek + 0.511e-3
!        call cdedxElec(e, -1, dedt)
!        write(*, *)sngl(ek), sngl(dedt)
!      enddo
!      end
     
!     ****************************************************************
!     *                                                              *
!     * cdedxElec:  gives -de/dx  (gev/(g/cm2) of e+/e-
!     *        if knock-on process is not to be included,            *
!     *        -de/dx by sternheimer is computed else                *
!     *        -de/dx by sternheimer - (-de/dx(recoil k.e > w) )     *
!     *        is computed.                                          *
!     *                                                              *
!     ************************ tested 87.09.19 ********************k.k
!
! /usage/  call cdedxElec(e, ic, dedt)

! -- input--
!        e: total energy of electron in GeV
!       ic: charge. -1 if electron. 1 if positron
!
!  *** note ***
!      Before calling this, cdedxEleci must have been called.
!
!
!
      subroutine cdedxElec(eini, ic, dedt)
       implicit none
!
!
#include  "ZdedxAir.h"
!
      real*8 eini, dedt
      integer ic
!      external cblkdedx

      real*8 emass, emassg
      real*8 ek, ein, e, gi, Beta2, x, a, b, c, cb, x0, x1
      real*8 dltx, wm, wlg
!
      real*8 peake/2.3d-7/, peak/0.858d0/
      data emass/0.511d0/, emassg/0.511d-3/


!
       if(jdef .eq. 0) then
          write(*, *) ' cdedxEleci must be called'
          stop 9999
       endif

       ek=eini-emassg
       ein=eini

       if( ek .lt. peake) then
              dedt=peak*sqrt(ek/peake)
       else
!                  energy in mev unit
              e=ein*1000.d0
              gi=emass/e
              Beta2= 1.d0 - gi**2
!                  x=log10(p/mc)
              x=log10(e**2 - 0.2611d0)/2.0d0 + 0.29d0
              a=stha
              b=sthb
              c=sthc
              cb=-c
              x0=sthx0
              x1=sthx1
              if(x .lt. x0) then
!                    4.605x - dlt
                  dltx=4.605d0*x
              elseif(x .lt. x1) then
                  dltx=cb - (x1-x)**3 * sthsa
              else
                  dltx=cb
              endif
              wm=e- emass
!              if(Knckon) then
!                Now we don't use Knckon, since big wm can do
!               the equivalent
                   if(wm .gt. w0) then
                       wlg=wlg0
                   else
                       wlg=log(wm)
                   endif
!              else
!                   wlg=log(wm)
!              endif
              if(ic .eq. -1) then
                  dedt=a/Beta2 * (b + 1.12d0 + wlg-Beta2 +dltx)
              else
                  dedt=a/Beta2 *( b+0.693d0 + wlg -2*Beta2 + dltx)
              endif
!                            convert it to gev/(g/cm^2)
              dedt=dedt *1.d-3
          endif
        end
