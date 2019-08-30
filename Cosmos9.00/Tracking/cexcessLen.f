!  cexcesLenght:
!         ************************************************************
!         *  compute exccess length by scattering of a charged ptcl
!         *  as compared to the straight pass
!         *
!                 !
!                 ! *
!                 ! *           dt given here is the excess path
!                 !   *         of * as compared with sqrt(r**2+l**2)-l
!                 !    *        it should be > sqrt(r**2+l**2)-l
!              l  !     *   
!                 !     *      i.e.  path length of * is sqrt(r**2+l**2) + dt
!                 !      *
!                     r
!               To get total time needed to run the * path, 
!             simply add dt/beta to the time in MovedTrack.
!    
      subroutine cexcessLen(dx, dy, dt)
      use modSetIntInf
      use modEMcontrol
      implicit none
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
#include  "ZmediaLoft.h"
      
      real*8  dx, dy  ! small dispalcement by scattering and geomag effect.
!          this treatment valid for only scattering displacement, but
!          for small displacement, we may add geomag. effect too.
! 
      real*8  dt   ! output in m  (not devided by beta)
      real*8 tmp, ee1, rr2, rp2, t, tang, tang1, za, za2
      real*8  dtp, en, u, g, al, tt, reall2, ge, reall
!
!           
      real*8 dtp0/.12/    !  dtp at rp2=0
      real*8  cthick2den
      real(8):: X0
!
      if(dx .eq. 0. .and.  dy .eq. 0.) then
         dt = 0.
      else
          za = TrackBefMove%pos%depth
          za2 = MovedTrack%pos%depth 
          ee1 = TrackBefMove%p%fm%p(4)
          reall = IntInfArray(ProcessNo)%length
          X0 = Media(MediaNo)%X0kg   ! r.l in Kg/m2
          t = IntInfArray(ProcessNo)%thickness/X0  ! r.l travelled 

          tmp=(Es/ee1)**2*cthick2den(.6666*za+ .3333*za2  )/X0
!
          rr2 = dx**2 + dy**2
          rp2 = rr2/tmp/reall**3
!                get coefficient for straight line for rp2 vs dtp
!               (average relation)
!                get tangent

          if(ee1 .gt. 500.e-3) then
             if(t .gt. 0.3) then
                tang=(0.49/(za+.1) + .053)*t + 1.13
             else
                tang=0.33 *t + 1.05
             endif
          else
             tang1=4.37e-2/ee1 + 1.3e-2
             tang=tang1*t+ 1.05
          endif
!
!                  get <dtp>
          dtp=tang*rp2+dtp0
!               get n for
!               dtp(>tt)/dtp(all)= tt**n/( .7544+tt**n)
!               this dtp is normalized by <dtp>
          en=19.375*rp2 +3.531
!               dtp= t**en/(beta+t**en) is a first good  approximation
!               which is solved for uniform random # as
!               t=( beta*u/(1-1))**(1/en). at small u,
!               t is under estimated so that it's correction is needed.
!               sample dtp ( <> normailzed )
          call rndc(u)
          if(u .lt. .3) then
!              g= (1.35-1.15)/6.5 * (en-10.)+ 1.35 =
             g= 3.077e-2*en +1.042
             al= ((3.33- 3.33/g)*u+1./g)/en
          else
             al=1./en
          endif

          tt=( 0.7544*u/(1.-u) )** al
!           convert to dt (in m)
          reall2=reall**2
          dt=dtp*tmp*reall2/2
!             it should be > sqrt(reall**2 + r**2) -reall
!             it may happen not so due to some approximation
          ge = sqrt(reall2 + rr2) -reall
          dt = max(ge, dt)
       endif
      end
