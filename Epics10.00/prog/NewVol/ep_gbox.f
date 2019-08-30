!
!        gbox:   general box.
!       Z                           
!       |   C     B     z-x side-1
!       |     
!       |         
!       O---------A-x
!
!
!            
!        
!      |   C'        B'  z-x side-2
!      |    
!      |      O'     A'  
!      |
!      -----------------
!
!             one vertex is on the origin.
!
!      The points A, B, C, A', B', C', O'
!        A, A', O' are on the xy plane.
!         A is on the X-axis.
!     each 4 vertex square must be coplaner within some accuracy
!
!   Data format in config is:
!       ox oy oz  Ax,  Bx, By, Bz, Cx, Cy, Cz, A'x, A'y,  
!                      B'x, B'y, B'z, C'x, C'y, C'z, O'x O'y 
!
!      where (ox,oy,oz) is the origin in the world coord.
!      
      subroutine eprgbox(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!         interface to read configuration data for "gbox"
!
       type(Component)::  comp  ! output. to recieve the config data.
 
       type(epPos)::  o, a, b, c, ap, bp, cp, op
       integer iax,  ibx, iby, ibz, icx, icy, icz, 
     *         iapx, iapy, ibpx, ibpy, ibpz, icpx, icpy, icpz,
     *         iopx, iopy
       parameter(iax=1, ibx=2, iby=3, ibz=4, icx=5, 
     *           icy=6, icz=7, iapx=8, iapy=9, ibpx=10, 
     *           ibpy=11, ibpz=12, icpx=13, icpy=14, icpz=15,
     *           iopx=16, iopy=17)


       logical flag

       flag = .false.

!           read gbox data as 'new-*'
!           gbox has 21 volume attributes and the direction cosines
!           of the 'x' and 'y'==> (1-6)
!
!             next is mandatory
        call eprpst(comp, 17, 17, 1, 6)
!
!           check some values
        a%x = Volat( comp%vol + iax)
        a%y = 0.
        a%z = 0.

        b%x = Volat( comp%vol + ibx)
        b%y = Volat( comp%vol + iby)
        b%z = Volat( comp%vol + ibz)

        c%x = Volat( comp%vol + icx)
        c%y = Volat( comp%vol + icy)
        c%z = Volat( comp%vol + icz)

        ap%x = Volat( comp%vol + iapx)
        ap%y = Volat( comp%vol + iapy)
        ap%z = 0. 

        bp%x = Volat( comp%vol + ibpx)
        bp%y = Volat( comp%vol + ibpy)
        bp%z = Volat( comp%vol + ibpz)

        cp%x = Volat( comp%vol + icpx)
        cp%y = Volat( comp%vol + icpy)
        cp%z = Volat( comp%vol + icpz)

        op%x = Volat( comp%vol + iopx)
        op%y = Volat( comp%vol + iopy)
        op%z = 0.
        o%x = 0.
        o%y = 0.
        o%z = 0. 
!                 z-x
        call epchkplane(o, a, b, c, 'OABC', flag)
!                 z-x 2
        call epchkplane(op, ap, bp, cp, "O'A'B'C'", flag)
!                 bottom
        call epchkplane(o, a, ap, op, "O'AA'O'", flag)
!             top
        call epchkplane(c, b, bp, cp, "CBB'C'", flag)
!              z-y  right
        call epchkplane(a, ap, bp, b, "AA'B'B", flag)
!             z-y near left
        call epchkplane(o, c, cp, op, "OCC'O'", flag)
        if(flag) then
           stop 9999
        endif
       end
!   ***************************************
      subroutine epchkplane(a, b, c, d, msg, flag)
      implicit none
#include "ZepPos.h"
!          check if 4 points a, b, c, d forms a plane
       type(epPos)::  a, b, c, d ! input. 4 points
      character*(*)  msg        ! input. message to be shown if
                                ! a,b,c,d are not coplaner
      logical flag              ! output.  made to be t, if
                                !       message was issued.
!            
       type(epPos)::  n1, n2, n3
      real*8 k1, k2, absn1, absn2, eps

      data eps/1.d-4/

      call ep3p2plane(a,  b, c, n1, k1)
      call epscalerProd(n1, n1, absn1)
      call ep3p2plane(a,  c, d, n2, k2)
      call epscalerProd(n2, n2, absn2)
      if(absn1 .eq. 0. .and. absn2 .eq. 0.) then
         call cerrorMsg(msg, 1)
         call cerrorMsg("4 points cannot form a plane", 1)
         flag = .true.
      elseif(absn1 .ne. 0.  .and. absn2 .ne. 0.) then
         call epvectorProd(n1, n2, n3)
         call epscalerProd(n3, n3, k1)
         if( k1 .gt. eps) then
            call cerrorMsg(msg, 1)
            call cerrorMsg("4 points cannot form a plane", 1)
            flag = .true.
         endif
      endif
      end         
      subroutine epbgbox(comp, pos, dir, length, icon)
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
                               !  by Volat( comp.vol + 1), etc
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real*8  length !  output length cm from pos to the boundary
       integer icon  ! output 0: length obtained. pos    is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
 
       real*8   l, xpa(2)
 
       type(epPos)::  o, a, b, c, ap, bp, cp, op
       integer iax,  ibx, iby, ibz, icx, icy, icz, 
     *         iapx, iapy, ibpx, ibpy, ibpz, icpx, icpy, icpz,
     *         iopx, iopy
       parameter(iax=1, ibx=2, iby=3, ibz=4, icx=5, 
     *           icy=6, icz=7, iapx=8, iapy=9, ibpx=10, 
     *           ibpy=11, ibpz=12, icpx=13, icpy=14, icpz=15,
     *           iopx=16, iopy=17)



       integer  jcon, np


        a%x = Volat( comp%vol + iax)
        a%y = 0.
        a%z = 0.

        b%x = Volat( comp%vol + ibx)
        b%y = Volat( comp%vol + iby)
        b%z = Volat( comp%vol + ibz)

        c%x = Volat( comp%vol + icx)
        c%y = Volat( comp%vol + icy)
        c%z = Volat( comp%vol + icz)

        ap%x = Volat( comp%vol + iapx)
        ap%y = Volat( comp%vol + iapy)
        ap%z = 0.

        bp%x = Volat( comp%vol + ibpx)
        bp%y = Volat( comp%vol + ibpy)
        bp%z = Volat( comp%vol + ibpz)

        cp%x = Volat( comp%vol + icpx)
        cp%y = Volat( comp%vol + icpy)
        cp%z = Volat( comp%vol + icpz)

        op%x = Volat( comp%vol + iopx)
        op%y = Volat( comp%vol + iopy)
        op%z = 0.
        o%x = 0.
        o%y = 0.
        o%z = 0. 

       
       np = 0
!              bottom 
       call epxpLand4vp(a, ap, op, o, pos, dir, l, icon, jcon)
       if(icon .le. 4 .and. l > 0.) then
          np = np +1
          xpa(np) = l
       endif

!            top
       call epxpLand4vp(b, bp, cp, c, pos, dir, l, icon, jcon)
       if(icon .le. 4 .and. l > 0.) then
          np = np +1
          xpa(np) = l
       endif
       if(np .eq. 2) goto 100

!           side x-z 1       
       call epxpLand4vp(o, a, b, c, pos, dir, l, icon, jcon)
       if(icon .le. 4 .and. l > 0.) then
          np = np +1
          xpa(np) = l
       endif
       if(np .eq. 2) goto 100

!           side x-z 2
       call epxpLand4vp(ap, bp, cp, op, pos, dir, l, icon, jcon)
       if(icon .le. 4 .and. l> 0.) then
          np = np +1
          xpa(np) = l
       endif
       if(np .eq. 2) goto 100

!           y-z 1       
       call epxpLand4vp(o, c, cp, op, pos, dir, l, icon, jcon)
       if(icon .le. 4 .and. l > 0.) then
          np = np +1
          xpa(np) = l
       endif
       if(np .eq. 2) goto 100
!           y-z 2
       call epxpLand4vp(a, ap, bp, b, pos, dir, l, icon, jcon)
       if(icon .le. 4 .and. l > 0. ) then
          np = np +1
          xpa(np) = l
       endif
       if(np .eq. 2) goto 100
!           not cross with the volume
       icon = -1
       goto 200
 100   continue
       if(xpa(1) .ge. 0. .and.  xpa(2) .ge. 0.) then
!             outside
          icon = 1
          length = min(xpa(1), xpa(2))
       elseif(xpa(1) .lt. 0. .and.  xpa(2) .lt. 0.) then
!             outside
          icon = -1
          length = max(xpa(1), xpa(2))
       else
!           inside
          icon = 0 
          length = max(xpa(1), xpa(2))
       endif
 200   continue
       end

!      **********************************
      subroutine epsgbox(comp, pos, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
!
!           judge if a given 'pos' is inside 'comp'
!         
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

 
       type(epPos)::  o, a, b, c, ap, bp, cp, op
       integer iax,  ibx, iby, ibz, icx, icy, icz, 
     *         iapx, iapy, ibpx, ibpy, ibpz, icpx, icpy, icpz,
     *         iopx, iopy
       parameter(iax=1, ibx=2, iby=3, ibz=4, icx=5, 
     *           icy=6, icz=7, iapx=8, iapy=9, ibpx=10, 
     *           ibpy=11, ibpz=12, icpx=13, icpy=14, icpz=15,
     *           iopx=16, iopy=17)



       type(epPos)::  dir
       real*8 length


        a%x = Volat( comp%vol + iax)
        a%y = 0.
        a%z = 0.

        b%x = Volat( comp%vol + ibx)
        b%y = Volat( comp%vol + iby)
        b%z = Volat( comp%vol + ibz)

        c%x = Volat( comp%vol + icx)
        c%y = Volat( comp%vol + icy)
        c%z = Volat( comp%vol + icz)

        ap%x = Volat( comp%vol + iapx)
        ap%y = Volat( comp%vol + iapy)
        ap%z = 0.

        bp%x = Volat( comp%vol + ibpx)
        bp%y = Volat( comp%vol + ibpy)
        bp%z = Volat( comp%vol + ibpz)

        cp%x = Volat( comp%vol + icpx)
        cp%y = Volat( comp%vol + icpy)
        cp%z = Volat( comp%vol + icpz)

        op%x = Volat( comp%vol + iopx)
        op%y = Volat( comp%vol + iopy)
        op%z = 0.
        o%x = 0.
        o%y = 0.
        o%z = 0. 

      if( pos%z .lt. 0.d0) then
         icon = 1
      elseif( pos%z .gt. max(b%z, c%z, bp%z, cp%z)) then
         icon = 1
      elseif(pos%x .lt. min(0.d0, c%x, cp%x, op%x) ) then
         icon = 1
      elseif(pos%x .gt.  max(a%x, ap%x, b%x, bp%x)) then
         icon = 1
      elseif(pos%y .lt. min(0.d0, a%y, b%y, c%y)) then
         icon = 1
      elseif(pos%y .gt. max(ap%y, bp%y, cp%y, op%y) ) then
         icon = 1
      else
!            draw half-line   with dir. (0,0,1)     
         dir%x = 0.
         dir%y = 0.
         dir%z = 1.d0
         call epbgbox(comp, pos, dir, length, icon)
         if(icon .ne. 0) then
            icon = 1
         endif
      endif

      end
!     **************************************
      subroutine epenvlpgbox(comp, org, abc)
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

 
       type(epPos)::  o, a, b, c, ap, bp, cp, op
       integer iax,  ibx, iby, ibz, icx, icy, icz, 
     *         iapx, iapy, ibpx, ibpy, ibpz, icpx, icpy, icpz,
     *         iopx, iopy
       parameter(iax=1, ibx=2, iby=3, ibz=4, icx=5, 
     *           icy=6, icz=7, iapx=8, iapy=9, ibpx=10, 
     *           ibpy=11, ibpz=12, icpx=13, icpy=14, icpz=15,
     *           iopx=16, iopy=17)






        a%x = Volat( comp%vol + iax)
        a%y = 0.
        a%z = 0.

        b%x = Volat( comp%vol + ibx)
        b%y = Volat( comp%vol + iby)
        b%z = Volat( comp%vol + ibz)

        c%x = Volat( comp%vol + icx)
        c%y = Volat( comp%vol + icy)
        c%z = Volat( comp%vol + icz)

        ap%x = Volat( comp%vol + iapx)
        ap%y = Volat( comp%vol + iapy)
        ap%z = 0.

        bp%x = Volat( comp%vol + ibpx)
        bp%y = Volat( comp%vol + ibpy)
        bp%z = Volat( comp%vol + ibpz)

        cp%x = Volat( comp%vol + icpx)
        cp%y = Volat( comp%vol + icpy)
        cp%z = Volat( comp%vol + icpz)

        op%x = Volat( comp%vol + iopx)
        op%y = Volat( comp%vol + iopy)
        op%z = 0.
        o%x = 0.
        o%y = 0.
        o%z = 0. 


      org%x = min(0.d0, c%x, cp%x, op%x)
      org%y = min(0.d0, a%y, b%y, c%y)
      org%z = 0.
      abc%x = max(a%x,  ap%x, bp%x, b%x) - org%x
      abc%y = max(ap%y,  bp%y, cp%y, op%y) - org%y
      abc%z = max(b%z, bp%z, cp%z, cp%z) - org%z

      NVTX = 8

      VTXx(1) = 0.
      VTXy(1) = 0.
      VTXz(1) = 0.
      
      VTXx(2) = a%x
      VTXy(2) = 0.
      VTXz(2) = 0.

      VTXx(3) = b%x
      VTXy(3) = b%y
      VTXz(3) = 0.

      VTXx(4) = c%x
      VTXy(4) = c%y
      VTXz(4) = 0.

      VTXx(5) = o%x
      VTXy(5) = o%y
      VTXz(5) = o%z

      VTXx(6) = ap%x
      VTXy(6) = ap%y
      VTXz(6) = ap%z

      VTXx(7) = bp%x
      VTXy(7) = bp%y
      VTXz(7) = bp%z

      VTXx(8) = cp%x
      VTXy(8) = cp%y
      VTXz(8) = cp%z

      end
!    *************************************
      subroutine epatlocgbox(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(17)
 
      integer i

      do i = 1, 17
         loc(i) = i
      enddo
      end


