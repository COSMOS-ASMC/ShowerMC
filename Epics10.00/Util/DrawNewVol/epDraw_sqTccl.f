      subroutine epDraw_sqTccl(comp, p, n)
      use sqTccl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               ! sqTccl in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  
      integer i
      real(8):: temp, temp2

      call epDraw_sqTccl0(comp, p, n)
      if( comp%struc == 'sqTccl' .or.
     *    comp%struc == 'sqTccl_w') then
      else
         do i = 1, n
            if( p(i)%x /= gpsep ) then
               call epc2v_sqTccl(comp, p(i), p(i))
            endif
         enddo
      endif
      end

      subroutine epDraw_sqTccl0(comp, p, n)
      use sqTccl
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               ! sqTccl in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  


      real(8):: fai, Dfai, temp
      real(8):: epsqTcclDfai
       type(epPos)::  p1
       type(epPos)::  p2, p3,p4


      real(8)::x
      real(8)::t1, t2
      integer n1,n2
      integer q
      logical kdgtest, drawopen
      call ep_sqTcclCnst(comp)


!         ccl
!      call epdrawCcl(r, h, thetamax, thetamin, p, n1)
!      n = n1
!      n = n + 1
!      p(n) = epPos(gpsep, 0.d0, 0.d0)
!      n = n + 1
!      p(n) = epPos(gpsep, 0.d0, 0.d0)
!      call epdrawCylEdg(p, n1, h,  p(n+1), n2)
!      n = n + n2

      fai = 0.
      n = 0
!//
      n1 = 0
!////////
      do while(.true.)
         Dfai =   epsqTcclDfai(fai)
         call epsqTcclGet4p(fai, Dfai, p1, p2, p3, p4, q)
         n =  n +1
         p(n) = p1
            n1 = n1 + 1            !!!
         n = n + 1
            n1 = n1 + 1            !!!
         p(n) = p1
         n = n + 1
            n1 = n1 + 1  !!!
         p(n) = p4
         fai = fai +Dfai
         if(fai > twopi-eps) exit
      enddo
!/////////
!      write(0,*) ' n1=',n1
!////////
      n = n + 1
      p(n) =epPos(gpsep, 0.d0, 0.d0)

      fai = 0.
!///////
      n2 = 0
!//////////
      do while(.true.)
         Dfai =   epsqTcclDfai(fai)
         call epsqTcclGet4p(fai, Dfai, p1, p2, p3, p4, q)
            n2 = n2 +1             !!!!
         n = n + 1
         p(n) = p2
         n = n + 1
            n2 = n2 +1  !!!!
         p(n) = p3
         n = n + 1
            n2 = n2 +1  !!!!
         p(n) = p3
         fai = fai + Dfai
         if(fai > twopi-eps) exit
      enddo
!//////
!      write(0,*) ' n2 =', n2
!////////
      n = n + 1
      p(n) = epPos(gpsep, 0.d0, 0.d0)
      n = n + 1
      p(n) = epPos(gpsep, 0.d0, 0.d0)

!!        bottom
      n = n + 1
      p(n) = epPos(-a, -b, 0.d0)
      n = n + 1
      p(n) = epPos(a,  -b, 0.d0)

      n = n + 1
      p(n) = epPos(gpsep, 0.0, 0.d0)


      n = n + 1
      p(n) = epPos(-a, b, 0.d0)
      n = n + 1
      p(n) = epPos(a,  b, 0.d0)
      
      n = n + 1
      p(n) = epPos(gpsep, 0.0, 0.d0)
      n = n + 1
      p(n) = epPos(gpsep, 0.0, 0.d0)

!         ///////////////////// cap of the ccl part
      t1 = thetamax*Torad
      t2 = thetamin*Torad
      fai = t1
      n1= 0
      n2 = n + 1

      call epsqTcclStepFai(fai, fai, Dfai)

      drawopen = abs( thetamin - thetamax ) .lt. 359.5
      do while(.true.)
         Dfai =   epsqTcclDfai(fai)
         if( fai > t2 ) exit
         call epsqTcclGet4p(fai, Dfai, p1, p2, p3, p4, q)
         n =  n +1
         p(n) = p1
            n1 = n1 + 1            !!!
         n = n + 1
            n1 = n1 + 1            !!!
         p(n) = p1
         n = n + 1
            n1 = n1 + 1  !!!
         p(n) = p4
         fai = fai +Dfai
      enddo
      n = n + 1
      p(n)%x = gpsep
!         draw  opened region if exists
      if( drawopen ) then
         p(n)%x = 0.
         p(n)%y = 0.
         p(n)%z = h
         n = n + 1
         p(n)= p(n2)
         n = n + 1
         p(n)%x = gpsep
      endif
      n = n + 1
      p(n)%x = gpsep

      call epdrawCylEdg(p(n2), n1, h,  p(n+1), n2)
      n = n + n2
      end
