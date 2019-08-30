!     ****************************************************************
!     *  This should not be used at kinetic energies < Mass          *
!     *                                                              *
!     * cdedxatHE: gives -de/dx (GeV/(kg/m**2)) of mu, pi, k, p
!     *        in air  (de/dx includes knockon electron energy)
!     *    Density effect depending on the air density is considered *
!     *
!     *        -de/dx by sternheimer is computed
!     *   ********** use prededx  for making const *******
!     *                                                              *
!     ************************ tested 87.09.19 ********************k.k
!
! /usage/  call cdedxInAir(ptcl, rhoin, dedt)
! -- input--
!     ptcl:  type ptcl
!    rhoin: density of air in kg/m**3
! -- output --
!     dedt; energy loss   GeV/(kg/m**2) 
!
!
!
         subroutine cdedxatHE(aPtcl, rhoin, dedt)
         implicit none
#include  "Zptcl.h"
#include  "ZdedxAir.h"
         type(ptcl):: aPtcl
         real*8  rhoin, dedt
!
!

         real*8 emass, emass2
         parameter (emass=.511d0, emass2=emass**2)
!
         real*8  ek, amass, a, b,  x0, c, x1, sa, rho, rl, e, gi
         real*8  pbmc, cb, dltx, p2, wm, wlg, x
!
!          see      sternheimer's consts. by p.r.b vol.3 (1971)3681
!                 to mev unit
       if(jdef .eq. 0) then
          write(*, *) ' cdedxEleci must be called'
          stop 9999
       endif

         if(aPtcl%charge .ne. 0) then
              ek=aPtcl%fm%p(4) *1.d3
              amass=aPtcl%mass * 1.d3
              a=7.68d-2
              b=17.87d0
              if(ek .lt. 50.d0*amass) then
                  x0=1.884d0
                  c=-10.97d0
                  x1=4.d0
                  sa=.247d0
              else
                  rho=min( max(rhoin, 1.d-7), 5.d-3)
                  rl=log10(rho)
                  x0= (((-.0165d0*rl-0.305d0)*rl-1.94d0)*rl
     *             -5.41d0)*rl-3.87d0
                  c=-4.0635d0+2.303d0*rl
                  if(rho .gt. 2.78605d-4) then
                      x1=4.d0
                      sa=((-.6872064d-01*rl-.5340530d0 )*rl
     *               -1.521159d0 )*rl-1.365158d0
                  else
                      x1=5.d0
                      sa=(((.04256d0*rl+0.7888d0)*rl+5.465d0)*rl
     *                   + 16.642d0)*rl+    18.855d0
                  endif
              endif
!                     energy in mev unit
              e=(ek+amass)
              gi=amass/e
              betasq = 1.d0 - gi**2
!                      x=log10(p/mc)
              pbmc =      (e/amass)**2 - 1.d0
              if(pbmc .gt. 0.) then
                  x=log10(pbmc)/2.0d0
                  cb=-c
                  if(x .lt. x0) then
!                        4.606x - dlt
                      dltx=4.606d0*x
                  elseif(x .lt. x1) then
                      dltx=cb - (x1-x)**3 * sa
                  else
                      dltx=cb
                  endif
                  p2=e**2 -amass**2
                  wm=2*emass*p2/( amass**2+ emass2+ emass*e*2)
!
                  if(wm .gt. w0) then
                     wlg = wlg0
                  else
                     wlg=log(wm)
                  endif
                  dedt=a/betasq *( b+0.693d0+wlg -2.d0*betasq + dltx)
!                                convert it to GeV/(kg/m2)
                  dedt=dedt *1.d-4 *aPtcl%charge**2
            else
                  dedt=1.d-9
            endif
         else
            dedt = 0.
         endif
      end
