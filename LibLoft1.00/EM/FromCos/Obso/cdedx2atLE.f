!#include "cblkdedx.h"
!      implicit none
!      real*8 etemp, amass, dedt, e
!      integer i
!
!      amass=    105.659e-3
!      amass=    938.659e-3
!      amass = 500.e-3
!      do  i=1, 200
!         e=      1.e-6* 10.**((i-1)/20.0)
!         etemp=e+amass
!         call cdedx2atLE(etemp,   amass, dedt)
!         write(*, *) sngl(e), sngl(dedt)
!      enddo
!      end
!     ****************************************************************
!     *                                                              *
!     * cdedx2atLE: gives -de/dx  (gev/(g/cm2)) of non electrons          *
!     *         this gives dedx at all energy region but
!     *         with constant air density
!     *                                                              *
!     ************************ tested 87.09.19 ********************k.k
!
! /usage/  call cdedx2atLE(e,amass, dedt)

! -- input--
!        e: total energy of ptcl     in GeV
!    amass: mass in GeV
! -- output --
!     dedt; energy loss / (g/cm2) in gev unit.
!
!
!
!
      subroutine cdedx2atLE(eini, amassi, dedt)
       implicit none
!
!
#include  "ZdedxAir.h"
      real*8 eini, amassi, dedt
      real*8  ek,  peak2, peak2e
!
      ek=eini-amassi
      if(ek .lt. 1.d-3) then    ! < 1MeV
         peak2e = 1.4d-6* (amassi/100.d-3)**1.6d0
         if(ek .lt. peak2e ) then
            if(ek .gt. 0.) then
               peak2 = 0.6d0*(ek/1.5d-5)**(-0.15d0)
               dedt=peak2*sqrt(ek/peak2e)
            else
               dedt = 0.
            endif
         else
            call cdedx2atLEa(eini, amassi, dedt)
         endif
      else
         call cdedx2atLEa(eini, amassi, dedt)
      endif
      if(dedt .lt. 0.) then
         dedt = 1.d-7
      endif
      end
!
      subroutine cdedx2atLEa(eini, amassi, dedt)
       implicit none
!
!
#include  "ZdedxAir.h"

!
      real*8 eini, amassi, dedt

      real*8 emass, emass2
      parameter(emass=.511d0, emass2=emass**2)
      real*8  e, amass, ein, gi, Beta2, x, a, b, c, cb, x0
      real*8 dltx, p2, wm, wlg,  x1


       if(jdef .eq. 0) then
          write(*, *) ' cdedxEleci must be called'
          stop 9999
       endif
       ein=eini
!
!                  energy in mev unit
          e=ein*1000.d0
          amass=amassi*1000.d0
          gi=amass/e
          if(gi .ge. 1.d0) then
              dedt=.01d0
          else
            Beta2= 1.d0 - gi**2
!                x=log10(p/mc)
            x=log10((e/amass)**2 - 1.d0) / 2.d0
            a=stha
            b=sthb
            c=sthc
            cb=-c
            x0=sthx0
            x1=sthx1
            if(x .lt. x0) then
!                  4.605x - dlt
                dltx=4.605d0*x
            elseif(x .lt. x1) then
                dltx=cb - (x1-x)**3 * sthsa
            else
                dltx=cb
            endif
            p2=e**2 -amass**2
            wm=2*emass*p2/( amass**2+ emass2+ emass*e*2)
            if(wm .gt. w0) then
               wlg = wlg0
            else
               wlg=log(wm)
            endif
            dedt=a/Beta2 *( b+0.693d0+wlg -2.d0*Beta2 + dltx)
        endif
        dedt=dedt *1.d-3
       end




