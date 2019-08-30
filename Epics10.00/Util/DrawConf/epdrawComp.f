      subroutine epdrawBox(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
!      integer how              ! input. i-th digit of how shows
                               !  whether the i-th box surface 
                               !  be drawn or not. i-th bit =0==>
                               !  not drawn, !=0==> drawn.
                               !  The surface is numbered as below.

!
!          |                                  1 x-y z=0
!          |     ******************           6 x-y z=c
!          |   * |          5    **           2 x-z y=0           
!          | *   |   6/        *  *           5 x-z y=b
!        c |********/ ********* 4 *           3 y-z x=0
!          | 3 b |/+++++++++++*+++*           4 y-z x=a 
!          |    /    1        *  *
!          |   /              * *
!          | /    2           *
!          |------------------------------
!                             a
!


       type(epPos)::  p(42)     ! output. (x,y,z) to describe
                               !   a box in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p
 

      logical kdgtest

      real*8 a, b, c


      n = 0
      a = Volat( comp%vol + boxa)
      b = Volat( comp%vol + boxb)
      c = Volat( comp%vol + boxc)


      if( kdgtest(how, 1) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
     
      if( kdgtest(how, 2) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = c

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = c

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 3) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = c

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = c

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 4) )then
         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = c

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = c

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 5) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = c

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = c

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 6) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = c

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = c

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = c

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = c

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      end
!     *********************
      subroutine epdrawPrism(comp,  p, n)
      use prism
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)
      integer,intent(in)::n

      integer i
       type(epPos)::  pp

      call epprismCnst(comp)
      call  epdrawPrism0(comp,  p, n)
      do i = 1, n
         if(p(i)%x /= gpsep ) then
            call epc2v_prism(comp, p(i), pp)
            p(i) = pp
         endif
      enddo
      end


      subroutine epdrawPrism0(comp,  p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
!      integer how              ! input. i-th digit of how shows
                               !  whether the i-th  surface 
                               !  be drawn or not. i-th bit =0==>
                               !  not drawn, !=0==> drawn.
                               !  The surface is numbered as below.

!
!          z                                  1 x-y z=0
!          |                                  2 x-z y=0
!          |                                  3 y-z slant
!          |                                  4 y-z slant
!          |   3     *& & & & &    5          5 x-z y=b 
!          |       *   *        *             
!          |     *        *   4   *
!          |   *    h 2      *      *
!          | *                  *     *b
!          |-----------------------*-------x
!                    c 1           a
!


       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !  a prism in local coordnate.
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p
 

      logical kdgtest

      real*8 a, b, c, h


      n = 0
      a = Volat( comp%vol + prisma)
      b = Volat( comp%vol + prismb)
      c = Volat( comp%vol + prismc)
      h = Volat( comp%vol + prismh)


      if( kdgtest(how, 1) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif
     
      if( kdgtest(how, 2) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = c
         p(n)%y = 0
         p(n)%z = h

         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 3) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = c
         p(n)%y = 0.
         p(n)%z = h

         n = n + 1
         p(n)%x = c
         p(n)%y = b
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 4) )then
         n = n + 1
         p(n)%x = a
         p(n)%y = 0.
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = c
         p(n)%y = 0.
         p(n)%z = h

         n = n + 1
         p(n)%x = c
         p(n)%y = b
         p(n)%z = h

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif


      if( kdgtest(how, 5) )then
         n = n + 1
         p(n)%x = 0.
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
!        -----------

         n = n + 1
         p(n)%x = c
         p(n)%y = b
         p(n)%z = h

         n = n + 1
         p(n)%x = a
         p(n)%y = b
         p(n)%z = 0.

         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      end
!     *****************************
      subroutine epdrawCyl(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a cyl in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  
      integer  i
       type(epPos)::  pp

      call epdrawCyl0(comp, p, n)

      if( comp%struc(1:5) == "cyl_y" .or. 
     *    comp%struc(1:5) == "cyl_x"  ) then
         do i = 1, n
            if( p(i)%x /= gpsep )   then
               call epc2v_cyl(comp, p(i), pp)
               p(i) = pp
            endif
         enddo
      endif
      end

      subroutine epdrawCyl0(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a cyl in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  


      logical kdgtest, drawopen

      integer n1, n2

      drawopen = abs( thetamin - thetamax ) .lt. 359.5
      call epdrawCcl(Volat( comp%vol + cylr), 0.d0,
!              min         max (not opposit !)
     *        thetamax, thetamin, p, n1)
!         draw  opened region if exists
      if( drawopen ) then
         p(n1)%x = 0.
         p(n1)%y = 0.
         p(n1)%z = 0.
         n1 = n1 + 1
         p(n1)= p(1)
         n1 = n1 + 1
         p(n1)%x = gpsep
      endif
      call epdrawCcl(Volat( comp%vol + cylr), Volat( comp%vol + cylh), 
     *        thetamax, thetamin, p(n1+1), n2)
      n = n1 + n2
      if( drawopen ) then
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z =  Volat( comp%vol + cylh)
         n = n + 1
         p(n)= p(n1+1)
         n = n + 1
         p(n)%x = gpsep
      endif
      n = n + 1
      p(n)%x = gpsep
!   
      if(kdgtest(howcyl, 1)) then
         call epdrawCylEdg(p, n1, 0.d0, p(n+1), n2)
         n = n + n2 
      endif
      if(kdgtest(howcyl, 2)) then
         call epdrawCylEdg( p(n1+1), n1, Volat( comp%vol + cylh),
     *        p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      end
      subroutine epdrawCylEdg(peri, n1, h,  p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
      integer n1    ! input. number of points in peri
       type(epPos)::  peri(n1)  ! input. the vertex of cyl-like object
                              !  (circle or part of it, etc      )
      real*8 h    !i nuput. place of the surface height.
                  !   h=0 for floor, h> 0 for ceil.
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
               ! a surface which cuts the cyl-like object.  The 
               ! surface makes the right angle with the axis
               ! of the cyl-like object at height h. (local coordinate).
               ! The size will be 2*n1+2 
      integer n ! output. number of vertex put in p.

       type(epPos)::  cent
      integer i

      cent%x = 0.
      cent%y = 0.
      cent%z = h

      do i = 1, n1-1
         p(i) = cent
      enddo
      n = n1 
      p(n)%x = gpsep

      do i = 1, n1
         p( n + i ) = peri(i)
      enddo
      n = n + n1
      n = n + 1
      p(n)%x = gpsep
      end
!     *****************************
      subroutine epdrawSphere(r0, zmin, zmax,  p, n)
      implicit none

#include "ZepDraw.h"
      real*8 r0               ! input.  radius of sphere
      real*8 zmin            ! input.  only z>=zmin  is drawn
      real*8 zmax            ! input.  only z<=zmax  is drawn 
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   a box in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension must be > 9*34
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      real*8 dt, theta, r, z
      integer i, n1, n1save, nvtxccl

      save n1save

      theta = asin(1.d0) * 2
      dt = theta / nvsph
      

      if(zmin .gt. -r0 .and.
     *     zmin .lt. r0  .and.
     *     zmin .lt. zmax ) then
         r = sqrt(r0**2 - zmin**2)
         call epdrawCcl(r, zmin, thetamax, thetamin, p, n)
         n1save = n
      else
         n = 0
      endif

      do i = 1, nvsph + 1
         z = r0 *  cos(theta)
         if(z .ge. zmin .and. z .le. zmax) then
!            if( i  .eq. 1  .or. i .eq. nvsph+1  ) then
!               not show 
!            else
               r = r0 *  sin(theta)
               call epdrawCcl(r, z, 
     *           thetamax, thetamin, p( n + 1 ), n1)
               n = n1 + n
               n1save = n1
!            endif
         endif
         theta = theta - dt
      enddo

      if(zmax .le.  r0 .and.
     *    zmax  .ge. -r0 .and.
     *    zmax .ge. zmin ) then
         r = sqrt(r0**2 - zmax**2)
         call epdrawCcl(r, zmax, thetamax, thetamin, p(n+1), n1)
         n1save = n1
      else
         n1 = 0
      endif
      n = n + n1
      

      n = n + 1
      p(n)%x = gpsep

      return
!     ***************
      entry epqsphere(nvtxccl)
!     ******************
      nvtxccl = n1save

      end
      subroutine epdrawPipe(comp, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   pipe in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      integer::i
       type(epPos)::  pp

      call epdrawPipe0(comp, p, n)
!     *******************************
      if( comp%struc(1:6) == "pipe_y" .or. 
     *    comp%struc(1:6) == "pipe_x"  ) then
         do i = 1, n
            if( p(i)%x /= gpsep )   then
               call epc2v_pipe(comp, p(i), pp)
               p(i) = pp
            endif
         enddo
      endif
      end

      subroutine epdrawPipe0(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   pipe in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      integer n1, n2
      integer lociru, locirl, locoru, locorl

      logical  kdgtest, drawopen
      
       type(epPos):: pls, pus, ple, pue

      drawopen = abs(thetamin - thetamax) .lt. 359.5

      call epdrawCcl(Volat( comp%vol + pipeir), 0.d0,
     *         thetamax, thetamin, p, n1)
      n = n1

      locirl = 1
      pls = p(1)       ! lower inner circle starting point
      ple = p(n1-1)    ! lower inner circle ending point

      call epdrawCcl(Volat( comp%vol + pipeir), 
     *  Volat( comp%vol + pipeh), thetamax, thetamin, p(n+1), n1)
      lociru = n+1
      pus = p(lociru)
      n = lociru + n1 -1
      pue = p(n-1)
      n = n + 1
      p(n)%x = gpsep
!     
!       outer wall
!
      if(drawopen) then
         n = n + 1
         p(n) = pls
      endif
      call epdrawCcl(Volat( comp%vol + pipeor), 0.d0,
     *      thetamax, thetamin,  p(n+1), n1)
      locorl = n+1
      n = n + n1
      if( drawopen ) then
         p(n) = ple
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n) = pus
      endif   
      call epdrawCcl(Volat( comp%vol + pipeor),
     *  Volat( comp%vol + pipeh), thetamax, thetamin,  p(n+1), n1)
      locoru = n+1
      n = n + n1 
      if( drawopen ) then
         p(n) = pue
         n = n + 1
         p(n)%x = gpsep
      endif   
      n = n +1
      p(n)%x = gpsep
!         draw the pipe cross-section part
      if(kdgtest(howcyl, 2)) then
         call epdrawPipeEdg(p(lociru), p(locoru), n1-1, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      if(kdgtest(howcyl,1) ) then
         call epdrawPipeEdg(p(locirl), p(locorl), n1-1, p(n+1), n2)
         n = n + n2
         n = n+ 1
         p(n)%x = gpsep
      endif
      end

      subroutine epdrawPipeEdg(p1, p2, n1, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

      integer n1     ! input. 
       type(epPos):: p1(n1), p2(n1)  ! input. ccl like vertex
       type(epPos):: p(*)  ! output.
      integer n      !  output. number of vertex put in p

      integer i

      do i = 1, n1
         p(i) = p1(i)
      enddo
      n =  n1 
      n = n+ 1
      p(n)%x = gpsep
      do i = 1, n1
         p(n+i) = p2(i)
      enddo
      n = n + n1
      n = n+ 1
      p(n)%x = gpsep
      end

!     *****************************
      subroutine epdrawCcl(r, h, t1, t2, p, n)
      implicit none

#include "ZepDraw.h"
#include "Zglobalc.h"
!
!          draw a circle of radius r.
!
      real*8 r                 ! input.  radius of the circle (cm)
      real*8 h                 ! input.  circle is at z=h. (cm)
      real*8 t1                ! input. azimuthal angle  region from t1 deg
      real*8 t2                ! input. to t2 deg  is drawn
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               !   a circle of radius r
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line
                               ! dimension must be as large as
                               ! nvccl+2 in ZepDraw.h
                               ! p(1) = p(last-1)
                               ! p(last) = gpsep (separator)
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      real*8 theta, dt, twopi, tmax
      integer i, imax

      
      twopi =  asin(1.d0)*4
      dt = twopi/nvccl 
      

!      theta = thetamax * Torad
      theta = t1 * Torad
!      imax = (thetamin - thetamax)/ 360.d0 * nvccl + 1
      imax = max(int( (t2-t1)/ 360.d0 * nvccl) ,1)
!      tmax = thetamin  * Torad
      tmax = t2  * Torad
!       adjust dt so that equi. dt
      dt = Torad*(t2-t1)/imax
      
!      i = 1
!      do while(.true.)
      n = 0
      do i = 1, imax+1
         n = n + 1
         p(i)%x = r*cos(theta)
         p(i)%y = r*sin(theta)
         p(i)%z = h
!         if(abs(theta-tmax) .le. dt/10.) goto 10
!         if(theta .eq. tmax) goto 10
         theta =min(theta + dt, tmax)
!         i = i + 1
      enddo
 10   continue
      n = n + 1
      p(n)%x = gpsep
      end
!       *******************
      subroutine epDraw_cap(comp, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   cap in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

       type(epPos)::  q1(nvccl+5), q2(nvccl+5)
      logical kdgtest
      integer ir, or, w1, h
      parameter( ir=1, or=2, w1=3, h=4)

      real*8 r, z

      integer n1, n2
      

      r = Volat( comp%vol + or)
      z = Volat( comp%vol + h)
      call epdrawSphere(r, z, r, p, n1)
      n = n1
      r = Volat( comp%vol + ir)
      call epdrawSphere(r, z, r, p(n+1), n1)
      n =  n + n1

      if(kdgtest(howcyl, 1)) then
         r =abs( Volat( comp%vol + w1))
         call epdrawCcl(r, z, thetamax, thetamin, q1, n1)
         r =
     *     sqrt( Volat( comp%vol + or)**2 - Volat( comp%vol + h)**2 )
         call epdrawCcl(r, z, thetamax, thetamin, q2, n1)
         call epdrawPipeEdg(q1, q2, n1-1, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      end

!       *******************
      subroutine epDraw_cone(comp, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   cone in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  


      integer ia, ib, ih, iap, ik
      parameter( ia = 1,  ib = 2,  ih = 3, iap=4, ik=5 )


      logical kdgtest
      integer n1, n2

      call epdrawElps(Volat( comp%vol + ia), 
     *   Volat( comp%vol + ib), 0.d0, thetamax, thetamin,   p, n1)
      n = n1
      call epdrawElps(Volat( comp%vol + iap), 
     *     Volat( comp%vol + ib)*Volat( comp%vol + ik),
     *     Volat( comp%vol + ih), thetamax, thetamin, p(n+1), n1)
      n = n + n1
      n = n + 1
      p(n)%x = gpsep

      if(kdgtest(howcyl, 1) ) then
         call epdrawCylEdg(p, n1, 0.d0, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif

      if(kdgtest(howcyl, 2) ) then
         call epdrawCylEdg(p(n1+1), n1, Volat( comp%vol + ih), 
     *               p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif

      end
!       *******************
      subroutine epDraw_ecyl(comp, p, n)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                               !   ecyl in local coordnate.

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  


      integer n1, n2
      logical kdgtest

      integer ra, rb, he
      parameter( ra=1, rb=2, he=3)

      real*8 a, b, h
!        
      a = Volat( comp%vol + ra)
      b = Volat( comp%vol + rb)
      h = Volat( comp%vol + he)

        
      call epdrawElps(a, b, 0.d0, thetamax, thetamin, p, n1)
      n = n1
      call epdrawElps(a, b, h, thetamax, thetamin, p(n+1), n1)
      n = n + n1
      n = n + 1
      p(n)%x = gpsep
!   
      if(kdgtest(howcyl, 1)) then
         call epdrawCylEdg(p, n1, 0.d0, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif

      if(kdgtest(howcyl, 2)) then
         call epdrawCylEdg(p(n1+1), n1, h, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif

      end
!     ********************
      subroutine epdrawElps(a, b, h, t1, t2, p, n)
      implicit none
#include "Zglobalc.h"

#include "ZepDraw.h"
      real*8 a    ! input. radius for x direction
      real*8 b    ! input. radius for y direction
      real*8 h    ! input. ellipse at hight h
      real*8 t1   ! input. azimuthal angle region from t1 deg to
      real*8 t2   ! input. t2 deg is drawn
       type(epPos)::  p(*)  ! output. vertex for ellipse
      integer n   ! output.  number of p.

      real*8  r
      integer i, nc
!          (x/a)^2 + (y/b)^2 = 1
!           tan A = y/x
!          r= sqrt( x^2 + y^2)  =|x|*sqrt( 1+ tan^2 A)
!         (x/a)^2 + (xtanA/b)^2 = 1
!         |x|= 1/sqrt ( (1/a)^2 + (tanA/b)^2)
      real*8 theta, dt, twopi, tmax, tant, minab
      twopi =  asin(1.d0)*4

      minab = min(a, b)
      dt = 360.d0/nvccl

      if( t2 < t1 ) then
         tmax = t2 + 360.d0
      else
         tmax = t2
      endif

      theta = t1
      n = 0
      do while(.true.)
         tant = tan(theta*Torad)
         r = sqrt(  (1.+ tant**2) / ( (1.0/a)**2 + (tant/b)**2) )
         n = n + 1
         p(n)%x = r*cos(theta*Torad)
         p(n)%y = r*sin(theta*Torad)
         p(n)%z = h
         if( theta == tmax ) exit
         theta =min( theta + dt*minab/r, tmax)
         if( theta > tmax - dt*minab/r/10.) then
            theta = tmax
         endif
      enddo
      n = n + 1
      p(n)%x = gpsep
      p(n)%y = gpsep
      p(n)%z = gpsep
      end
