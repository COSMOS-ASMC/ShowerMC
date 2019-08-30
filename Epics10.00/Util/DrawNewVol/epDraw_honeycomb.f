      module drawhoneycom

      integer,save::toomany_y=30    !if ny > toomany_y, next is applied for ydirecion
      integer,save::toomany_x=15    !if nx > toomany_x, next is applied for xdirecion
      real(8),save::drawportion=0.2  ! only 20 % of periheral region is drawn
                          ! if this is > 0.5, all honeycom will be drawn

      end module drawhoneycom

      subroutine epDraw_honeycomb(comp, p, n)
      use honeycomb
      use drawhoneycom
      implicit none
#include "ZepManager.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               ! honeycomb in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  

      real(8):: temp, temp2
      integer:: i
      logical,save:: first=.true.
      integer:: icon

      if(first) then
!  This is special for honeycomb: 
!     avoid center part drawing ; how much of the cener  part
!     is fixed by this  data.
         call copenf(iowk, "Data/honeycomb.dat", icon)
         if( icon ==0 ) then
            read(iowk,*) toomany_x, toomany_y, drawportion
            close(iowk)
         else
            write(0,*) 'Warning: honeycomb.dat cannot be opened'
            write(0,*) 'default values are employed'
            toomany_x=15
            toomany_y=30
            drawportion=0.15
         endif
         first =.false.
      endif

      call epDraw_honeyc0(comp, p, n)

!      if( comp.struc == "honeycomb" ) then
!      elseif( comp.struc == "honeycomb_x" ) then
!      elseif( comp.struc == "honeycomb_xy" ) then
!      elseif( comp.struc == "honeycomb_yz" ) then
!      elseif( comp.struc == "honeycomb_y" ) then
!      elseif( comp.struc == "honeycomb_z" ) then
         do i = 1, n
           if( p(i)%x /= gpsep ) then
              call epc2v_honeycomb(comp, p(i), p(i))
!              temp = b- p(i).y
!              temp2 = p(i).x
!              p(i) = epPos( temp, temp2, p(i).z)
           endif
        enddo
!      elseif( comp.struc == "honeycomb_zx" ) then
!         do i = 1, n
!           if( p(i).x /= gpsep ) then
!              call epc2v_honeycomb(comp, p(i), p(i))
!              temp = p(i).x
!              temp2 = b-p(i).y
!              p(i) = epPos( p(i).z, temp2, temp)
!           endif
!        enddo
!      else
!         write(0,*) 'comp.struc=', comp.struc, ' invalid'
!         write(0,*) 'detected in epDraw_honeycomb'
!         stop
!      endif
      end subroutine  epDraw_honeycomb
   
      subroutine epDraw_honeyc0(comp, p, n)
      use honeycomb
      use drawhoneycom
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)      ! output. (x,y,z) to describe
                               ! honeycom in local coordnate.
                               ! (x,y,z)= gpsep is a separator
                               ! to be converted to a blank line

      integer  n               ! output.  number of (x,y,z) data
                               ! put in p.  
      integer nx, ny
      integer nx1, nx2, ny1, ny2
      integer k, i, j
      real(8)::xp(4)  ! nodal point of x
      real(8)::x, y, z
      real(8)::array(2)  ! array(2) will by y for given x=array(1)

      call ep_honeycmCnst(comp)

      nx =  (Xs + a)/Sx
      ny =  (Ys + b)/HSy
!      write(0,*) ' Sx =', Sx, ' HSy=', HSy
!      write(0,*) ' nx =', nx, ' ny=', ny
      nx1 = nx
      nx2 = nx1+1
      ny1 = ny
      ny2 = ny1+1
      if( nx >= toomany_x )  then
         if(drawportion < 0.5 ) then
            nx1 = nx*drawportion
            nx2 = nx*(1.-drawportion)
         endif
      endif
      if(ny >= toomany_y) then
         ny1 = ny*drawportion
         ny2 = ny*(1.-drawportion)
      endif
      n = 0
      do i = 0,  ny
!           at z=0
         do j = 0, nx
            if( i > ny1 .and.   i < ny2 .and.
     *          j > nx1 .and.   j< nx2 ) cycle
            xp(1) = -HLS + j*Sx
            xp(2) = xp(1) + LS
            xp(3) = xp(2) +  LIcost
            xp(4) = xp(3) + LS
            do k = 1, 4
               array = nodal_y(xp(k), i)
               if(array(1) < Xs .or. array(1) > Xs+a ) cycle
               if(array(2) < Ys .or. array(2) > Ys+b ) cycle 
               n = n + 1
               p(n) = epPos(array(1), array(2), 0.d0)
            enddo
         enddo

         n = n + 1
         p(n) = epPos(gpsep, 0.d0, 0.d0)
!           at z= h
         do j = 0, nx
            if( i > ny1 .and.   i < ny2 .and.
     *          j > nx1 .and.   j< nx2 ) cycle
            xp(1) = -HLS + j*Sx
            xp(2) = xp(1) + LS
            xp(3) = xp(2) +  LIcost
            xp(4) = xp(3) + LS
            do k = 1, 4
               array = nodal_y(xp(k), i)
               if(array(1) < Xs .or. array(1) > Xs+a ) cycle
               if(array(2) < Ys .or. array(2) > Ys+b ) cycle 
               n = n + 1
               p(n) = epPos(array(1), array(2), h)
            enddo
         enddo
         n = n + 1
         p(n) = epPos(gpsep, 0.d0, 0.d0)
         n = n + 1
         p(n) = epPos( gpsep, 0.d0, 0.d0)
      enddo
      end   subroutine epDraw_honeyc0
