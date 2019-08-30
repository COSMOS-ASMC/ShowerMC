      subroutine epDraw_cpipe(comp, p, n)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   a cpipe in local coordnate.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line
                              ! dimension of p must be >+ (nvccl+2)*2
      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      integer ih,  isa, iea, iir, ior
      integer ix1, iy1, ix2, iy2
      parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
      parameter (ix1=5, iy1=6, ix2=7, iy2=8)


      real*8   sa, ea


      real*8 t1, t2, sai, eai,  sa1, sa2, ea1, ea2
      integer  n2,  case
      logical isinside
      logical isinside2
      real*8 x
      isinside(x) = mod(thetamin-thetamax+360.d0, 360.d0) .ge.
     *               mod(x-thetamax+360.d0, 360.d0)
      isinside2(x) = mod(ea-sa+360.d0, 360.d0) .ge.
     *               mod(x-sa+360.d0, 360.d0)



      call epangcpipe(comp, 2, sa, ea)
      call epangcpipe(comp, 1, sai, eai)

      if(sai .eq. 0. .and. eai .eq. 0.) then
         case = 1
         sa1 = sa
         ea1 = ea
         sa2 = sa1
         ea2 = ea1
      elseif(sai .eq. 0. .and. eai .gt. 359.9999d0 
     *        .and. eai .lt. 360.d0) then
         case = 4
         sa1 =mod( (sa+ea)/2+360.0d0, 360.d0)
         ea1 = sa1+360.d0
         sa2 = sa
         ea2 = ea
         if(sa2 .lt. sa1) sa2 = sa2 + 360.d0
         if(ea2 .lt. sa2) ea2 = ea2 + 360.d0
      elseif( isinside2(sai) ) then
         case = 2
         sa1 = sa
         ea1 = ea
         sa2 = sai
         ea2 = eai
         if(sa2 .lt. sa1) sa2 = sa2 + 360.d0
         if(ea2 .lt. sa2) ea2 = ea2 + 360.d0
      else
         case = 3
         sa1 = sai
         ea1 = eai
         sa2 = sa
         ea2 = ea
         if(sa2 .lt. sa1) sa2 = sa2 + 360.d0
         if(ea2 .lt. sa2) ea2 = ea2 + 360.d0
      endif

      n = 0
      if(.not. isinside(sa1) .and. isinside(ea1) ) then
         t1 = thetamax
         t2 = ea1
      elseif(isinside(sa1) .and. isinside(ea1)) then
         if(isinside2(thetamin)) then
!              there is two regions
            t1 =  sa1
            t2 = mod(thetamin,360.d0)
            if(t2 .lt. t1) t2 = t2+360.d0
            call epcpipewall(comp, case, t1, t2, sa2, ea2, p(n+1), n2)
            n = n+ n2
            t1 =mod(thetamax, 360.d0)
            t2 = ea1
            if(t2 .lt. t1) t2 = t2+360.d0
            call epcpipewall(comp, case, t1, t2, sa2, ea2,
     *           p(n+1), n2)
            n = n + n2
            goto 100
         else
            t1 = sa1
            t2 = ea1
         endif
      elseif(isinside2(thetamax) .and. isinside2(thetamin)) then
         t1 =mod(thetamax, 360.d0)
         t2 =mod(thetamin, 360.d0)
      elseif(isinside2(thetamin) .and. isinside(sa1) ) then
         t1 = sa1
         t2 = mod(thetamin, 360.d0)
      else
         return                 !  no part to be drawn
      endif

      if(t2 .lt. t1) t2 = t2 + 360.d0
      call epcpipewall(comp, case, t1, t2, sa2, ea2,
     *     p(n+1), n2)
      n = n + n2
 100  continue
      end

      subroutine epcpipewall(comp, case, t1, t2, sa2, ea2,  p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component):: comp
      integer case  !  input 1~4
      real*8 t1, t2  ! input. starting and ending angle
      real*8 sa2, ea2  ! input. starting and ending angles for nodal points

       type(epPos)::  p(*)    !  output. points for the wall .



      integer  n             ! output.  number of (x,y,z) data

      integer ih,  isa, iea, iir, ior
      integer ix1, iy1, ix2, iy2
      parameter (ih = 2,  isa=3, iea=4, iir=9, ior=1)
      parameter (ix1=5, iy1=6, ix2=7, iy2=8)



       type(epPos)::  pos, dir
      integer loc1ru, loc1rl, loc2ru, loc2rl, nsv, nv, i, j
      real*8   dt, eps, temp, r, theta, tmax, length
      integer nc1, nc2, icon, n2
      logical kdgtest

      data eps/1.d-5/


      dt = 360.d0/nvccl
      theta = t1
      tmax =  t2
      if(tmax .lt. t1) then
         tmax= t2 + 360.
      endif
         
      if(case .le. 2) then
         r = Volat( comp%vol + ior)
      else
         r = Volat( comp%vol + iir)
      endif
      nc1 = 0
      nc2 = 0
      i = 0
      do while(.true.)
         if( theta .ge. sa2  .and. nc1 .eq. 0) then
            theta = sa2
            nc1 = nc1 +1
         elseif( theta+dt/4. .gt. sa2 .and. nc1  .eq. 0  .and. 
     *           i .gt. 0 ) then
            theta = sa2
            nc1 = nc1 +1
         endif

         i = i + 1
         p(i)%x = r*cos(theta*Torad)
         p(i)%y = r*sin(theta*Torad)
         p(i)%z = 0.d0
         if(theta .eq. tmax) goto 10
         theta =min(theta + dt, tmax)
         if(theta .gt. ea2 .and. nc2 .eq. 0) then
            theta = ea2
            nc2 = nc2 + 1
         elseif( theta+dt/4. .gt. ea2 .and. nc2  .eq. 0) then
            theta = ea2
            nc2 = nc2 +1
         endif
      enddo
 10   continue
      nv = i
      n = i + 1
      p(n)%x = gpsep
      loc1rl = 1

      i = n
      nc1 = 0
      nc2 = 0
      theta = t1
      loc1ru = n+1
      do j = 1, nv
         p( n +j)%x = p(j)%x
         p( n +j)%y = p(j)%y
         p(n+j)%z = Volat( comp%vol + ih)
      enddo
      n = n + nv
      n = n+1
      p(n)%x = gpsep

      n = n + 1
      p(n)%x = gpsep
!     ------------   
      
      nsv = n 
      loc2rl = n + 1

      do i = 1, nv
         dir = p(i)
         dir%z = 0.
!                   normalize
         temp =sqrt( dir%x**2 + dir%y**2 )
         if(temp .gt.  0.) then
            dir%x = dir%x/temp
            dir%y = dir%y/temp
            dir%z = 0.
         else
            dir%x = 1.
            dir%y = 0.
            dir%z = 0.
         endif

         if(case .le. 2) then
!               from outer circle to inner circle
            dir%x= -dir%x
            dir%y= -dir%y
         endif
         pos%x = p(i)%x + eps *dir%x
         pos%y = p(i)%y + eps *dir%y
         pos%z = eps
         call epbcpipe(comp, pos, dir, length, icon)
         if(icon .eq. -1) then
            p(n+i)%x = p(i)%x
            p(n+i)%y = p(i)%y
            p(n+i)%z = p(i)%z
         else
            p(n+i)%x = pos%x + length*dir%x
            p(n+i)%y = pos%y + length*dir%y
            p(n+i)%z = 0.
         endif
      enddo
      n = n + nv
      n = n+1
      p(n)%x = gpsep

      loc2ru = n+1
      do i = 1, nv
         p(n+i) = p(nsv+i)
         p(n+i)%z = Volat( comp%vol + ih)
      enddo
      n = n+ nv
      n = n+1
      p(n)%x = gpsep
      n = n+1
      p(n)%x = gpsep
      
      
      if(kdgtest(howcyl, 1)) then
!         call epdrawPipeEdg(p(loc1ru), p(loc2ru), nv+1, p(n+1), n2)
         call epdrawPipeEdg(p(loc1ru), p(loc2ru), nv, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      if(kdgtest(howcyl, 2)) then
!         call epdrawPipeEdg(p(loc1rl), p(loc2rl), nv+1, p(n+1), n2)
         call epdrawPipeEdg(p(loc1rl), p(loc2rl), nv, p(n+1), n2)
         n = n + n2
         n = n + 1
         p(n)%x = gpsep
      endif
      end

