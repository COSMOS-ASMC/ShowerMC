!         Compute geometrical relation in the spherical env.
!
!                |
!                |
!                |               /               *          .
!                |              /              * t    .
!                |             /             *.  
!                |            /         .  * A
!                |           / t'   .    *
!                |          /   . L    *         zenith angle t, t'
!                |      B  /.        *           at vertical height h and h'
!                |      . /        *  H=R0+h     Radius of the earth R0
!             ^  |  ^.   / H'    *               length AB =L
!        ^       |      /  ^   *            
!                |     /     *^       H cos(t) - H'cos(t')=L  
!                | R0 /    *    ^     H sin(t) = H'sin(t')    
!                |   /   *        ^   H' = sqrt(L**2 +H**2 - 2LHcos(t))
!                |  /  *          ^        
!                | / *                   
!                |/                 ^  cos(t')= (Hcos(t) - L)/H'
!                                        
!                                   ^         
!                                      sin(t')=Hsin(t)/H'
!                                        
!
!        test program.
!
!      program testutil
!      implicit none
!      include '../../Zglobalc.h'
!      include 'Zearth.h'
!      
!      real*8  H, sh, cosz, L, hp, cnewcos, cnewsin, cnewh, cost, sint
!      
!      do sh = 0., 30.d3, 10000.
!         H = Eradius + sh
!         do cosz=1., 0., -0.2
!            do L=0., 100.d3, 1000.
!               cost = cnewcos(H, cosz, L)
!               sint = cnewsin(H, cosz, L)
!               if(abs(cost**2 + sint**2 -1.d0) .gt. 1.d-7) then
!                  write(*,*) ' ***** '
!               endif 
!               write(*, *) cost**2 + sint**2 
!               hp = cnewh(H, cosz, L) - Eradius
!               write(*, *)' h=', sh, ' cosz=', cosz, ' L=', L, ' hp=',hp
!            enddo
!         enddo
!      enddo
!        another test; for clenbetween2h
!+++++++++++++++++++++++++++++++++++++++++++++++++++
!      program testutil
!      implicit none
!#include "Zearth.h"
!         real*8 clenbetween2h, h1, h2, cost, ans, cnewcos
!         real*8 costp
!         h1 = Eradius + 0.d3
!         h2 = Eradius + 100.d3
!         cost = 0.7
!         ans = clenbetween2h(h1, h2, cost)
!         write(*,*) ans
!         costp = cnewcos(h1, cost, ans)
!         ans = clenbetween2h(h2, h1, -costp)
!         write(*,*) ans
!       end
!c 
!     ----------------get cos(t')------------------------    
!        see also cnewcossin 
      real*8 function cnewcos(H, cost, L)
! 
      implicit none
      real*8 H, cost, L
!
      real*8 eps/1.d-8/,  tmp
      
      tmp = L/H
      if(tmp .lt. eps) then
          cnewcos = (cost - tmp) / 
     *      ( ( tmp * (1.-cost**2)/2 -cost)*tmp +1.)
      elseif(abs(cost) .ne. 1.d0) then
            cnewcos = (cost - tmp)/ sqrt( (tmp - cost*2)*tmp +1.)
      else
            cnewcos = cost
      endif
      end
!     ---------------get sin(t') -----------------------
      real*8 function cnewsin(H, cost, L)
      implicit none
      real*8 H, cost, L
!
      real*8 cnewh,  sint
      
      sint = sqrt(1.d0 - cost**2)
      cnewsin = H * sint/ cnewh(H, cost, L)
      end
!     ------------- get H'----------------------------
      real*8 function cnewh(H, cost, L)
      implicit none
      real*8 H, cost, L, tmp

      real*8 eps/1.d-8/

      tmp = L/H
      if(tmp .lt. eps) then
         cnewh = H * ( (tmp * (1.d0 -cost**2)/2 - cost )* tmp + 1.d0)
      else
         cnewh =H* sqrt( ( tmp - cost*2)*tmp + 1.d0 )
      endif
      end
!     -----------------------------------
      subroutine cnewcossin(h1, cos1, leng, h2, cos2, sin2)
      implicit none
      real*8  h1  !  input.  radial distance from the earth center
      real*8  cos1  ! input.  cos of zenith angle at h1
      real*8 leng    ! input.  length in m along cos1 from h1
      real*8  h2    ! output.  radial  distance from the earth center
      real*8  cos2  ! outpu.  cos of zenith angle at h2.
      real*8  sin2  ! output.  sin of //

      real*8 sin1, cnewh

!        h1 cos1 - h2 cos2 = leng
!        h1 sin1 = h2 sin2
      
      sin1 = sqrt(1.d0 - cos1**2)
      h2 = cnewh(h1, cos1, leng)
      sin2 = h1 * sin1/h2 
      cos2 = (h1*cos1 -leng)/h2
      end
!
!     ------------- get L----------------------------
!       get length between two points with radius h1
!       and h2. The zenith angle at h1 is cost.
!       the length has a sing, so that we can reach
!       h2 by proceeding the given length 
!       from h1 with the zenith angle cost
       real*8 function clenbetween2h(h1, h2, cost)
       implicit none
       real*8 h1, h2, cost

       real*8 sint, costp, sintp
       character*120 text

       sint = sqrt(1.d0 - cost**2)
       sintp = h1*sint /h2
       if(sintp .le. 1.0d0) then

          costp = sqrt(1.d0 - sintp**2)*sign(1.d0, cost)
       else
          if(abs(1.d0-sintp**2) .lt. 1.d-6) then
             costp = 0.
          else
             write(text, *) 'h1, h2, cost=', h1, h2, cost, 
     *          ' sintp=',sintp
             call cerrorMsg(text, 1)
             call
     *       cerrorMsg('h1,h2,cost invalid at clenbetwee2h', 0)
          endif 
       endif
!
       clenbetween2h = h1* cost - h2 * costp
      end
!
!     ------------- get L--- subroutine version with
!                            return condition ---------
!       get length between two points with radius h1
!       and h2. The zenith angle at h1 is cost.
!       the length has a sing, so that we can reach
!       h2 by proceeding the given length 
!       from h1 with the zenith angle cost
!       In some cases, the given condition cannot be
!       satisfied. In that case, icon is set to be 1
!       Such a case would happen if the track is far away
!       from the Earth center.
       subroutine clenbetw2h(h1, h2, cost, leng, icon)
       implicit none
       real*8 h1, h2, cost
       real*8 leng ! output
       integer icon ! output

       real*8 sint, costp, sintp
       character*120 text

       sint = sqrt(1.d0 - cost**2)
       sintp = h1*sint /h2
       icon = 0
       if(sintp .le. 1.0d0) then
          costp = sqrt(1.d0 - sintp**2)*sign(1.d0, cost)
          leng = h1* cost - h2 * costp
       else
          if(abs(1.d0-sintp**2) .lt. 1.d-6) then
             costp = 0.
             leng = h1* cost 
          else
             icon = 1
          endif 
       endif
      end
