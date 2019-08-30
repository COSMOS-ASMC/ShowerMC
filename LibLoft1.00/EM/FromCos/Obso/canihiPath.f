!     ****************************************************************
!     *
!     * canihiPath:  e+ annihilaion prob. per r.l
!     * capanihiea:  samples energy of 2 gammas and angles
!     *
!
!
!  **** note ***
!       before calling canihiea, canihip must be called and
!

!
      subroutine canihiPath(ein, prob,  path)
      implicit none
#include "Zmass.h"
      real*8 ein  !   input. positron energy in Gev
      real*8 prob !   output. prob. per r%l
      real*8 path !   output. sampled path in r%l
!

      real*8 g, csc, g2m, srg2m
!
      real*8  u
!            constm=.3*z/a*x0ing = media.basearea*2
!            constm=.03*z/a*x0inkgpm2
!
!            in  $elmag must have been fixed beforehand.

       real*8 constm/5.475/
       real*8 geqv1/1.00000001d0/
       save geqv1, constm
!
      g=ein/masele
      if( g .ge. geqv1) then
         if(g .gt. 2.5 .and. g .lt. 25.) then
!             approx within 1%
            csc=(((((0.9382535e-07*g-.8791180e-05)*g+0.3338919e-03)*g
     *          -.6609205e-02)*g
     *          +0.7308900e-01)*g-.4516664)*g+ 1.534200
         elseif(g .lt. 2.5) then
            g2m=g**2-1.
            srg2m=sqrt(g2m)
            csc=( ((g+4.)*g+1.)/g2m * log(g+srg2m) - (g+3.)/srg2m)
     *           /(g+1.)
         else
            csc=1.6* g**(-7./9.)
         endif
!         prob = csc * media.basearea
         prob = csc * constm
         call rndc(u)
         path = - log(u)/prob
      else
         prob = 1.d10
         path = 0.
      endif
      end
!
!      ************
      subroutine canihiea(ein, eg1, eg2, cos1, cos2)
!      ************
      implicit none
#include "Zmass.h"
      real*8 ein
      real*8 eg1  !  produced gamma energy of higher one.
      real*8 eg2  !                            lower one.
      real*8 cos1 ! cos angle of eg1
      real*8 cos2 ! cos angle of eg2

      real*8 g,  gc, bc,  sints, sint2
!
      real*8  u,  egs,  x, x1, x2, a,  aq, c2, a22
      real*8 geqv1/1.00000001d0/
      integer count
      save geqv1

      count=0
       g=ein/masele
       if(g .le.  geqv1) then
          x = 0.5d0
       else
          a = g + 1.0
          c2 = a + 2.7182*g/a
          x1 =1./(a + sqrt(g**2 -1.))
          x2 = 0.5d0
          a22 = a**2-2.
          aq = a22 + 2*a
          do while(.true.)
             call rndc(u)
             x = x1 * exp(log( (1.0d0-x1)/x1 )*u)
             call rndc(u)
             count = count +1
             if( count .gt. 100 ) then
                write(0,*) '****** loop in canihiea'
                write(0,*) '****** ein=',ein
                write(0,*) ' ein=Me forced '
                x=0.5
                ein = masele
                g = 1.
                goto 100
             endif
             if(u .le.  (aq -a**2*x -1./x)/a22) goto 100
          enddo
 100      continue
       endif
       eg2 = x * ( ein + masele )
       eg1 = ein -  eg2 + masele
!
       eg1 = max(eg1, eg2)
       eg2 = ein - eg1 + masele
       gc=sqrt( (g+1.)/2 )
       bc=sqrt((g-1.)/(g+1.))
       egs=masele*gc
       if(g .gt. geqv1) then
!          cms cos of gammma 1
          u  =max(min( (eg1/gc/egs -1.0)/bc, 1.0d0), 0.d0)
       else
          call rndc(u)
       endif
!          cms sin^2
       sints=1.d0-u**2
!          lab sin^2 of gamma1
       sint2= (egs/eg1)**2 *sints
       if(sint2 .gt. 1.d0) then
          cos1 = 0.
       else
          cos1 = sqrt(1.d0 - sint2)
       endif
!          lab sin^2 of gamma2
       sint2= (egs/eg2)**2 * sints
       if( sint2 .gt. 1.d0) then
          cos2 = 0.
       else
          cos2=sqrt(1.d0 - sint2)
       endif
!
!           pcos = gcE*(bc + cos*)
!           cos*= -u so cos< 0 if bc+cos* < 0  
!
       if(bc-u .lt. 0.) then
         cos2=-cos2
       endif
       end
