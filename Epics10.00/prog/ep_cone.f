!
!        cone:
!
!        (rather a megaphone with an elliptical cross-section)
!        If the ellipse of one side is made to be very small, this
!        becomes essentially a cone.
!        The canonical form assumes the zaxis be the center line.
!
!
!
!                            ^   ** 
!                   ^          * | *
!               ^             *  |   *
!           *^               *   a'   * 
!        a * *               *   |    *
!       Y  *o*- - -  - - - - - - -    *------> Z
!          *|*               *        *
!           * _               *      *
!           |       _          *    *
!           |               -     **
!           |
!           X
!     
!      a: longer radius of the ellipse at z=0
!      b: shorter radius of the ellipse at z=0
!     a': longer radius of the ellipse at z=h
!
!       (but actually  could be a <= b, a'<=a)
!   Data format in config is:
!       ox oy oz  a  b h  a'  x.x x.y x.z  y.x y.y y.z
!
!      where (ox,oy,oz) is the origin in the world coord.
!            (a b h a') is the parameters to describe the cone
!            (x.x ..  y.z) optional. direction cosines of the
!      'a'-axis and 'b'-axis in the world coordinate. If not
!       given, (1,0,0) and (0,1,0) are assumed.
!      
      subroutine eprcone(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
!
!         interface to read configuration data for "cone"
!
       type(Component)::  comp  ! output. to recieve the config data.
       character*120 msg
 
       integer ia, ib, ih, iap, ik
       parameter( ia = 1,  ib = 2,  ih = 3, iap=4, ik=5 )

       real*8 a, b, h, ap
!
!           read cone data as 'new-*'
!           cone has 4 volume attributes and the direction cosines
!           of the 'a' and 'b'==> (1-6)
!
!             next is mandatory
        call eprpst(comp, 4, 5, 1, 6)
!
!           next is optional
!           check some values
        a = Volat( comp%vol+ia )
        b = Volat( comp%vol +ib )
        h = Volat( comp%vol +ih )
        ap = Volat( comp%vol +iap )
        if(a  .le. 0. .or. b .le. 0. .or. h .le. 0. 
     *      .or. ap .le. 0.) then
           write(msg, *) comp%cn, '-th component: a=', a,
     *    ' b=', b, ' h=', h, " a'=",ap, 
     *    ' for cone;  invalid(must>0)'
           call cerrorMsg(msg, 0)
        endif
        Volat( comp%vol +ik ) = ap/a        ! embed k. 
       end
!   ***************************************
      subroutine epbcone(comp, pos, dir, length, icon)
       implicit none
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

       type(Component):: comp  ! input. you can extract volume parameters
                               !  by Volat( comp.vol +1), etc
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real*8  length !  output length cm from pos to the boundary
       integer icon  ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
 
       integer ia, ib, ih, iap, ik
       parameter( ia = 1,  ib = 2,  ih = 3, iap=4, ik=5 )

       real*8 a0, b0, h,  k


       type(epPos)::  cp
       integer nb, jcon, kcon
!
       integer i
!
       real*8 p, q, r, alfa, d, leng(4), temp, eps

#ifdef DEBUG
       debug = -1
#endif

       nb = 0
       a0 = Volat( comp%vol +ia)
       b0 = Volat( comp%vol +ib)
       h = Volat( comp%vol +ih)
       k = Volat( comp%vol +ik)

       alfa = (k - 1.d0)/h
!
!       solve:   p *l^2 + 2q * l + r = 0,  where
!
       p = (b0 * dir%x)**2 +  (a0*dir%y)**2 -(a0*b0*alfa*dir%z)**2
       q = pos%x*dir%x*b0**2 + pos%y*dir%y*a0**2 -
     *      (a0*b0)**2 * alfa*dir%z*( 1.d0 + alfa* pos%z)
!     *      2*(a0*b0)**2 * alfa*dir.z*( 1.d0 + alfa* pos.z)

       r = (pos%x*b0)**2 + (pos%y*a0)**2 - (a0*b0 *
     *      (1.d0 +  alfa* pos%z ))**2


!
!
       d = q**2 - p*r
       if(p .eq. 0.d0 .and.  q .ne. 0.) then
          length = -r/q/2
          if(length .ge. 0.) then
             nb = nb + 1
             leng(nb) = length
          endif
       elseif(p .ne. 0.d0  .and. d .ge. 0.d0) then
!             p != 0  d>=0
          temp = sqrt(d)
          length = (-q + temp)/p
          if(length  .ge. 0.) then
             nb = nb + 1
             leng(nb) = length
          endif
          length = (-q - temp)/p
          if(length  .ge. 0.) then
             nb = nb + 1
             leng(nb) = length
          endif

       endif
       if(dir%z  .ne. 0.) then
          length = (h - pos%z)/dir%z
          if(length  .ge. 0.) then
             nb = nb + 1
             leng(nb) = length
          endif
          length = -pos%z/dir%z
          if(length  .ge. 0.) then
             nb = nb + 1
             leng(nb) = length
          endif
       endif
       if( nb  .eq. 0) then
          icon = -1
       else
          call epscone(comp, pos, jcon)
!
!            eps is to assure that the boundary is
!            inside the volume.
!
          if( jcon .eq. 0 ) then
             eps =  -EpsLeng
          else
             eps =  EpsLeng
          endif
          length = 1.d30
          do i = 1, nb
             cp%x = (leng(i)+eps)*dir%x +  pos%x
             cp%y = (leng(i)+eps)*dir%y +  pos%y
             cp%z = (leng(i)+eps)*dir%z +  pos%z
             call epscone(comp, cp, kcon)
             if(kcon .eq. 0) then
                if(leng(i) .lt. length) then
                   length = leng(i)
#ifdef DEBUG
                   debug = i
#endif
                endif
             endif
          enddo
          if(length .eq. 1.d30) then
             icon = -1
          else
             icon = jcon
          endif
       endif
       end
!      **********************************
      subroutine epscone(comp, pos, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         

       integer ia, ib, ih, iap, ik
       parameter( ia = 1,  ib = 2,  ih = 3, iap=4, ik=5 )

       real*8 a0, b0, h, k
       
       real*8 a, b,  alfa


       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside



      h = Volat( comp%vol +ih)

      if( pos%z .lt. 0.d0 ) then
         icon = 1
      elseif( pos%z .gt. h ) then
         icon = 1
      else
         k = Volat( comp%vol +ik)
         alfa = (k-1.d0)/h
         a0 =  Volat( comp%vol +ia)
         b0 =  Volat( comp%vol +ib)
         a  = a0*(1.d0 + alfa*pos%z)
         b  = b0*(1.d0 + alfa*pos%z)
         if(a .eq. 0. or. b .eq. 0.) then
            icon = 1
         elseif( (pos%x/a)**2 + (pos%y/b)**2 .gt. 1.d0) then
            icon =1 
         else
            icon = 0
         endif
      endif
      end
!     **************************************
      subroutine epenvlpcone(comp, org, abc)
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

 
       integer ia, ib, ih, iap, ik
       parameter( ia = 1,  ib = 2,  ih = 3, iap=4, ik=5 )

       real*8 k
!
      
       
      k = Volat( comp%vol +ik) 
      org%x = min(-Volat( comp%vol +ia), -Volat( comp%vol +iap))
      org%y = min(-Volat( comp%vol +ib), -Volat( comp%vol +ib)*k)
      org%z = 0.d0
      abc%x = abs(org%x)*2
      abc%y = abs(org%y)*2
      abc%z = Volat( comp%vol +ih)

      end
!     *************************************
      subroutine epatloccone(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer  loc(*)
      integer i
      do i = 1, 5
         loc(i) = i
      enddo
      end

