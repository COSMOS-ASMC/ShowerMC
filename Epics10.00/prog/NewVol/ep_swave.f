!  swavey board or corrugated sheet.



!               +  / +
!   z        +    /    +
!   |     +      / d     +
!   |           /          +                       +
!   |        *  *           /+                   +
!   |     *   |    *       /   +-------------- +        *
!   |   *     | h    *    /d     x            x       *
!   | *       |        * /         x       x        *
!   |*        |         *            x   x      L *
!   |---------------------------------------------|
!   | |                   *                   *   |  wavy part
!   | |                    *                 *    |
!   | |                     x              x      |  base box
!   | |o                       x        x         |
!   | |                           x  x            |
!   | |                                           |
! 0 | ----|---------------------------------------|--
!   0     s                                       e
!
!   Data format in config is:
!       ox oy oz  L h o  d s e  [dir]
!
!      where (ox,oy,oz) is the origin in the world coord.
!         L, h,  d,  e > 0  s>=0
!         o >= h
!
!      
      subroutine eprswave(comp)
       implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "swave"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*150 msg

       integer iL, ih, io, id, is, ie 
       parameter( iL = 1,  ih = 2,  io=3, id=4, is= 5, ie=6)

       real*8 L, h, o, d, s, e
!
!           read swave data as 'new-*'
!           swave has 5 volume attributes and the direction cosines
!           (1~6)
!
!             next is mandatory
        call eprpst(comp, 6, 6, 1, 6)
!
        L = Volat( comp%vol + iL)
        h = Volat( comp%vol + ih)
        o = Volat( comp%vol + io)
        d = Volat( comp%vol + id)
        s = Volat( comp%vol + is)
        e = Volat( comp%vol + ie)

        if(L  .le. 0. .or. h .le. 0. .or. o .lt. 0.
     *    .or. d .le. 0. .or. e .le. 0. .or. s .lt. 0. 
     *    .or. (o .lt.  h) ) then
           write(msg, *)
     *      comp%cn, '-th component: some of L=', L,
     *    ' h=', h, ' o=',o, ' d=',d, ' s=',s,' e=',e,
     *    ' for swave;  invalid'
           call cerrorMsg(msg, 1)
           if(o .lt. h) then
              call cerrorMsg('o< h', 0)
           endif
        endif
       end
!   ***************************************
      subroutine epbswave(comp, pos, dir, length, icon)
      implicit none
#include "Zglobalc.h"
#include "ZepTrackp.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"

!
!        find length to the boundary of 'comp' from 'pos'
!        with direction cos 'dir'
!     'pos' and 'dir' are given in this 'comp' local coordinate.
! 

       type(Component):: comp    ! input. you can extract volume parameters
                               !  by Volat( comp.vol + 1), etc
       type(epPos)::  pos        ! input.  position.
       type(epDirec)::  dir      ! input. direction cosinse

      real*8  length            !  output length cm from pos to the boundary
      integer icon              ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume

      integer iL, ih, io, id, is, ie 
      parameter( iL = 1,  ih = 2,  io=3, id=4, is= 5, ie=6)

      real*8  el
 

      integer nx
      real*8 leng, lenga(4)


       type(epPos)::  org, b 
       type(ep3Vec)::  abc 

      integer jcon, n0, nb, nmin, nmax, kcon, i, j
      real*8 xc, yc, zc, temp
      real*8 ep_swavef
      external ep_swavef

      real*8 eps
      save eps
!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, o
      common /Zswave/  L, h, d, s, e, a, hL, qL, o
!     *************************************
      data eps/1.d-9/

!
!        get bounding box

      call epenvlpswave(comp, org, abc)
!         next has been set in above call. and set in /Zswave/
!      L = Volat( comp.vol + iL)
!      hL = L/2
!      qL = L/4
!      h = Volat( comp.vol + ih)
!      o = Volat( comp.vol + io)
!      d = Volat( comp.vol + id)
!      s = Volat( comp.vol + is)
!      e = Volat( comp.vol + ie)
!        
!      a = 16.0*h/L**2


!         obtain the crossing point with the bounding box
      call kxplbx( pos%x-org%x, pos%y-org%y, pos%z - org%z,
     *            dir%x, dir%y, dir%z, abc%x, abc%y, abc%z,
     *            el, jcon)
      if(jcon .eq. -1) then
         icon = -1
         return                 ! ********
      endif

      
!       position n0
      n0 = int( pos%x/hL )

!       crossing point with the bounding box
      nb = int( ( pos%x + dir%x*el )/hL )

      nmin = int((s+eps)/hL)
      nmax = int((e-eps)/hL)

!       see if pos is inside 
      if(jcon .eq. 0) then
!            already in bounding box. so see if in swave part.
         if( ep_swavef(comp,1, pos%x) .ge. pos%z .and.
!     *       ep_swavef(comp,0, pos.x) .le. pos.z ) then
     *       pos%z  .ge. 0.0d0 ) then
            icon = 0
         else
            icon = 1
         endif
      else
         icon = 1
      endif
!           counter for the crossing points
      nx = 0
      if(icon .eq. 0) then

!            inside; crossing point with the upper swave
         call epswavecross(n0, pos, dir, 1, leng, kcon)
         if(kcon .eq. 0) then
            nx = nx + 1
            lenga(nx) = leng
         else
            if(dir%x .gt. 0.) then
               j = 1
            else
               j = - 1 
            endif
            do i = n0+j, nb, j
               call epswavecross(i, pos, dir, 1,  leng, kcon)
               if(kcon .eq. 0) then
                  nx = nx +1
                  lenga(nx) = leng
                  goto 10
               endif
            enddo
 10         continue
         endif
!               inside; crossing point with the lower swave
         call epswavecross(n0, pos,  dir, 0,  leng, kcon)
         if(kcon .eq. 0) then
            nx = nx + 1
            lenga(nx) = leng
!         else
!            if(dir.x .gt. 0.) then
!               j = 1
!            else
!               j = - 1 
!            endif
!            do i = n0+j, nb, j
!               call epswavecross(i, pos, dir, 0,  leng, kcon)
!               if(kcon .eq. 0) then
!                  nx = nx +1
!                  lenga(nx) = leng
!                  goto 20
!               endif
!            enddo
 20         continue
         endif

!              examin front/back wall
         if( dir%y .lt. 0.) then
!                 front wall
            leng = -pos%y/dir%y

            xc = pos%x + dir%x * leng
            zc = pos%z + dir%z * leng
            if(xc .ge. s .and. xc .le. e) then

               temp = ep_swavef(comp,1,xc)

!               if( temp .ge. zc .and. temp-t .le. zc) then
               if( temp .ge. zc .and. zc .ge. 0.d0) then
                  nx = nx + 1
                  lenga(nx) = leng
               endif
            endif
         elseif(dir%y .gt. 0.) then
!               back wall
            leng = (d - pos%y)/dir%y
            xc = pos%x + dir%x * leng
            zc = pos%z + dir%z * leng

            if(xc .ge. s .and. xc .le. e) then
               temp = ep_swavef(comp,1,xc)
!               if( temp .ge. zc .and. temp-t .le. zc) then
               if( temp .ge. zc .and.  zc .ge. 0.d0) then
                  nx = nx + 1
                  lenga(nx) = leng
               endif
            endif
         endif
!               side wall
         if(dir%x .gt. 0.) then
            leng = ( e - pos%x )/dir%x
            yc = pos%y + dir%y * leng
            zc = pos%z + dir%z * leng
            if(yc .ge.  0.d0 .and. yc .le. d) then
               temp = ep_swavef(comp,1, e)
!               if( temp .ge. zc .and. temp-t .le. zc) then
               if( temp .ge. zc .and.  zc .ge. 0.d0) then
                  nx = nx + 1
                  lenga(nx) = leng
               endif
            endif
         elseif(dir%x .lt. 0.) then
            leng = ( s - pos%x )/dir%x
            yc = pos%y + dir%y * leng
            zc = pos%z + dir%z * leng
            if( yc .ge. 0. .and. yc .le. d) then
               temp = ep_swavef(comp,1, s)
!               if( temp .ge. zc .and. temp-t .le. zc) then
               if( temp .ge. zc .and.  zc .ge. 0.d0 ) then
                  nx = nx + 1
                  lenga(nx) = leng
               endif
            endif
         endif
!/////////
         if(nx .eq. 0) then
            write(0,*) '*********  error ; nocrossing'
            write(0,*) 'pos%x,y,z=',pos%x,pos%y,pos%z,
     *              ' dir%x,y,z=',dir%x,dir%y, dir%z
            write(0,*) ' org%x,y,z=',org%x,org%y,org%z,
     *       ' abc%x,y,z=',abc%x, abc%y, abc%z
            write(0,*) ' nb=',nb
            debug = 1

            call epswavecross(n0, pos,  dir, 0,  leng, kcon)
            call epswavecross(n0, pos,  dir, 1,  leng, kcon)
            call epswavecross(n0-1, pos,  dir, 0,  leng, kcon)
            write(0,*) ' n0=',int( pos%x/hL )
            write(0,*) ' pos%x=',pos%x, ' hL=',hL
            write(0,*) ' epU=',ep_swavef(comp,1, pos%x)
            write(0,*) ' epL=',ep_swavef(comp,0, pos%x)
            write(0,*) ' pos%z=',pos%z
            stop
         endif
      else
!          pos is outside
         if(jcon .eq. 1)  then 
!              outside of bounding box
!             get crossing point with the bounding box
            b%x = pos%x + dir%x * el
            b%y = pos%y + dir%y * el
            b%z = pos%z + dir%z * el
         else
            b%x = pos%x
            b%y = pos%y
            b%z = pos%z
            el = 0.
         endif

         n0 = int(b%x/hL)
         temp =ep_swavef(comp,1, b%x)
         if( temp .le. b%z ) then
!               xp should be with upper swave
            if(dir%x .ne. 0.)  then
               call epswavecross(n0, b, dir, 1, leng, kcon)
               if(kcon .eq. 0) then
                  nx = nx + 1
                  lenga(nx) = leng + el
               elseif(dir%x .gt. 0. .and. n0+1 .le. nmax) then
                  call epswavecross(n0+1, b, dir, 1, leng, kcon)
                  if(kcon .eq. 0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
               elseif(dir%x .lt. 0. .and. n0-1 .ge. nmin) then
                  call epswavecross(n0-1, b, dir, 1, leng, kcon)
                  if(kcon .eq. 0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
	       endif	
            endif
!         elseif( temp-t .ge. b.z ) then   
         elseif( b%z .ge. 0.d0 ) then   
!               xp should be with lower swave
            if(dir%x .ne. 0.)  then
               call epswavecross(n0, b, dir, 0, leng, kcon)
               if(kcon .eq. 0) then
                  nx = nx + 1
                  lenga(nx) = leng + el
               elseif(dir%x .gt. 0. .and. n0+1 .le. nmax) then
                  call epswavecross(n0+1, b, dir, 0, leng, kcon)
                  if(kcon .eq. 0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
               elseif(dir%x .lt. 0. .and. n0-1 .ge. nmin) then
                  call epswavecross(n0-1, b, dir, 1, leng, kcon)
                  if(kcon .eq. 0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
               endif
            endif
         endif

!             examine front and back wall

         if( dir%y .ne. 0. ) then
!                 front wall
            leng = - b%y/dir%y
            if(leng .ge. 0.) then
               xc = b%x + dir%x * leng
               zc = b%z + dir%z * leng
               if(xc .ge. s .and. xc .le. e) then
                  temp = ep_swavef(comp,1,xc)
!                  if( temp .ge. zc .and. temp-t .le. zc) then
                  if( temp .ge. zc .and.  zc .ge. 0.d0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
               endif
            endif
!               back wall
            leng = (d - b%y)/dir%y
            if(leng .ge. 0.) then
               xc = b%x + dir%x * leng
               zc = b%z + dir%z * leng

               if(xc .ge. s .and. xc .le. e) then
                  temp = ep_swavef(comp,1,xc)
!                  if( temp .ge. zc .and. temp-t .le. zc) then
                  if( temp .ge. zc .and. zc .ge. 0.d0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
               endif
            endif
         endif
!               side wall
         if(dir%x .ne. 0.) then
            leng = ( e - b%x )/dir%x
            if(leng .ge. 0.) then
               yc = b%y + dir%y * leng
               zc = b%z + dir%z * leng
               if(yc .ge.  0.d0 .and. yc .le. d) then
                  temp = ep_swavef(comp,1,e)
!                  if( temp .ge. zc .and. temp-t .le. zc) then
                  if( temp .ge. zc .and. zc .ge. 0.d0) then
                     nx = nx + 1
                     lenga(nx) = leng + el
                  endif
               endif
            endif
            leng = ( s - b%x )/dir%x
            if(leng .ge. 0.) then
               yc = b%y + dir%y * leng
               zc = b%z + dir%z * leng
               if( yc .ge. 0. .and. yc .le. d) then
                  temp = ep_swavef(comp,1,s)
!                  if( temp .ge. zc .and. temp-t .le. zc) then
                  if( temp .ge. zc .and. zc .ge. 0.d0) then
                     nx = nx + 1
                     lenga(nx) = leng +el
                  endif
               endif
            endif
         endif
      endif
      if(nx .eq. 0) then
         icon = -1
      else
         length = lenga(1)
         do i = 2, nx
            if(length .gt. lenga(i)) then
               length = lenga(i)
            endif
         enddo
      endif            
      end          
      subroutine epswavecross(n, pos, dir, ul, leng, kcon)
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "ZepPos.h"
#include "Zcnfig.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
!
      integer n   ! input. n-th half swave part is examined
      integer ul  !  0 --> lower swavey part 1--> upper wavy part.
       type(epPos)::  pos
       type(epDirec)::  dir

      real*8 leng
      integer kcon

      real*8 aa, bb, cc, dd, cx, cy

!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, o
      common /Zswave/  L, h, d, s, e, a, hL, qL, o
!     *************************************
      real*8 thick 


      if(ul .eq. 1) then
         thick = o  ! not 0
         if(dir%x .eq. 0. .and. dir%z .eq. 0.) then
            kcon = -1
         else
            aa = dir%x**2
            bb = (2*pos%x - (2*n+1)*hL)*dir%x  + (-1)**n*dir%z/a
            cc = (pos%x - hL*n)*(pos%x - (n+1)*hL) +
     *           (-1)**n * (pos%z-thick)/a
            if(aa .eq. 0.) then
               leng = -cc /bb
               if(leng .gt. 0.) then
                  kcon = 0
                  return        ! ******
               endif
            else
               dd = bb**2 - 4*aa*cc
               if(dd .lt. 0.) then
                  kcon = -1
                  return        ! *******
               else
                  dd = sqrt(dd)
                  leng = (-bb - dd)/2/aa
                  if(leng .gt. 0.) then
                     cx = pos%x + dir%x*leng 
                     cy = pos%y + dir%y*leng 
!                       cz = pos.z + dir.z*leng 
                     if( cx .ge. max(hL*n,s) .and.
     *                    cx .le. min(hL*(n+1),e)) then
                        if( cy .ge. 0. .and. cy .le. d) then
                           kcon = 0
                           return ! *******
                        endif
                     endif
                     kcon = -1
                  endif

                  leng = (-bb + dd)/2/aa
                  if(leng .gt. 0.) then
                     cx = pos%x + dir%x*leng 
                     cy = pos%y + dir%y*leng 
                     if( cx .ge. max(hL*n, s) .and. 
     *                    cx .le. min(hL*(n+1),e) ) then
                        if( cy .ge. 0. .and. cy .le. d) then
                           kcon = 0
                           return !  *********
                        endif
                     endif
                     kcon = -1
                  else
                     kcon = -1
                  endif
               endif
               kcon = -1
            endif
         endif
      else
!           bottom 
         if(dir%z .ne. 0.)  then
            leng = - pos%z/dir%z
            if(leng .ge. 0.d0) then
               cx = pos%x + leng*dir%x
               cy = pos%y + leng*dir%y
               if(cx .ge. s .and. cx .le. e) then
                  if(cy .ge. 0. .and. cy .le. d) then
                     kcon = 0
                     return !   *****
                  endif
               endif
            endif
            kcon = -1
         endif
         kcon = -1
      endif            
      end

      subroutine epsswave(comp, pos, icon)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "ZepDirec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

       type(epPos)::  org 
       type(ep3Vec)::  abc 

!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, o
      common /Zswave/  L, h, d, s, e, a, hL, qL, o
!     *************************************



      integer  n0
      real*8  ep_swavef


      external ep_swavef

      call epenvlpswave(comp, org, abc)
!           common has been set by the above call

      if( pos%x .lt. org%x .or. pos%x .gt. org%x + abc%x) then
         icon  = 1
      elseif(pos%y .lt. org%y .or. pos%y .gt. org%y + abc%y) then
         icon = 1
      elseif(pos%z .lt. org%z .or. pos%z .gt. org%z + abc%z) then
         icon = 1
      else
         n0 = int( pos%x/hL )

         if( ep_swavef(comp,1,pos%x) .ge. pos%z .and.
!     *       ep_swavef(comp,0,pos.x) .le. pos.z ) then
     *        pos%z .ge. 0.d0 ) then
            icon = 0
         else
            icon = 1
         endif
      endif
      end
!     **************************************
      subroutine epenvlpswave(comp, org, abc)
      implicit none

#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

!
!        give the envloping box of the component
!
       type(Component)::  comp  ! input.   component.
       type(epPos)::  org       ! output.  origin of the enveloping box
                               !          in local coord. 
       type(ep3Vec)::  abc      ! output.  a,b,c of the box

!         L, h, t, d,  e > 0  s>=0
      integer iL, ih, io, id, is, ie 
      parameter( iL = 1,  ih = 2,  io=3, id=4, is= 5, ie=6)


!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, o
      common /Zswave/  L, h, d, s, e, a, hL, qL, o
!     *************************************

      real*8  zmin, zmax
      integer n
      real*8  ep_swavef

      external ep_swavef
!
      call ep_swavesetcom(comp)

      n = int((s-qL)/L)
      if( s .le. (qL+n*L) .and. e .ge. (qL+n*L)) then
         zmax = h+o
      else
         zmax = max( ep_swavef(comp,1,s), ep_swavef(comp,1,e))
      endif
!      n = int( (s-(qL+hL))/L ) 
!      if( s .le. (qL+hL) + n*L .and. e .ge.(qL+hL)+n*L )then
!         zmin = -h-t
!      else
!         zmin = min( ep_swavef(comp,0, s),  ep_swavef(comp,0,e))
!      endif

      zmin = 0.d0

      org%x = s
      org%y = 0.d0
      org%z = zmin
      abc%x = (e-s)
      abc%y = d
      abc%z = zmax - zmin
      NVTX = 0
      end
!     *************************************
      subroutine epatlocswave(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(6)
 
      integer i

      do i = 1, 6
         loc(i) = i
      enddo
      end
!     _____________________________ 
      subroutine ep_swavesetcom(comp)
!      set consts in the common area
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


!
       type(Component)::  comp  ! input.   component.


!         L, h, t, d,  e > 0  s>=0
      integer iL, ih, io, id, is, ie 
      parameter( iL = 1,  ih = 2,  io=3, id=4, is= 5, ie=6)


!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, o
      common /Zswave/  L, h, d, s, e, a, hL, qL, o
!     *************************************

!
      L = Volat( comp%vol + iL)
      hL = L/2
      h = Volat( comp%vol + ih)
      o = Volat( comp%vol + io)
      d = Volat( comp%vol + id)
      s = Volat( comp%vol + is)
      e = Volat( comp%vol + ie)
      qL = L/4

      a = 16.0*h/L**2

      end

      real*8 function ep_swavef(comp, ul, x)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

!
!      compute swavey function values at x
!      this gives 0.0 vale if x is outside of the 
!      defined region
!       
       type(Component)::  comp  ! input.   component.
      integer ul  ! input. 0:  lower wavy part.
                  !        1:  upper wavy part.
      real*8   x  ! input. x position  
!
!
!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, o
      common /Zswave/  L, h, d, s, e, a, hL, qL, o
!     *************************************

      real*8 f, g
      integer n
      integer cnsave
      save cnsave
      data cnsave/-1/

      f(x) = -a*x*(x-hL)  ! basic funciton at x=0 to hL. (half swave leng).
!                 max is at qL. the value is h. 
      g(x,n) = (-1)**n * f(x-hL*n)  ! n is integer to specify which 
                                    ! swave (in unit of half swave). 0~hL is 0
                                    !     hL~ L is 1....

      if(comp%cn .ne. cnsave) then
         call ep_swavesetcom(comp)
         cnsave = comp%cn
      endif

      n =int( x/hL )


      if(ul .eq. 1) then
         ep_swavef = g(x, n)
         ep_swavef = ep_swavef +o
      else
         ep_swavef = 0.d0
      endif

      end
