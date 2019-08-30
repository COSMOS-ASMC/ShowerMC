      subroutine eprcap(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp  ! output. to recieve the config data.
       character*120 msg
 
       integer ir, or, w1, h
       parameter( ir=1, or=2, w1=3, h=4)

       real*8 r1, r2, w
!
!           read Cap configuration data as 'new-1'
!            cap has 3 volume  attributes, its z-axis direction cos
!            (7-9) are to be read.
        call eprpst(comp, 3, 4, 7, 9)

!           check some values
!((((((((((((
        r1 = Volat( comp%vol +ir )
        r2 = Volat( comp%vol +or )
        w  = Volat( comp%vol +w1 )
!))))))))))))
        if(r1 .ge. r2) then
           write(msg, *) comp%cn, '-th component: r1=', r1,
     *    ' >= r2 =', r2, ' for Cap'
           call cerrorMsg(msg, 0)
        endif
        if(abs(w) .gt. r1 ) then
!               |w1| > r1
           write(msg, *) comp%cn,'-th component: w1=',w,
     *    ' >= r1 =', r1, ' for Cap'

           call cerrorMsg(msg, 0)
        endif
!           compute h for later use and embed it in
!           vol(4) which is not used yet.
!           vol(1) = r1, vol(2)=r2, vol(3)=w1
!           vol(4)= h  
!                           *   *    
!                       *   .   .    *
!                    *   .         .    *
!                  *   .              .   *   
!                 <-----------x------->
!                 \     w2    | w1=w /
!                    \       h|    /
!                       \     |  / r1
!                     r2   \  |/
!                            \x origin
! 
!((((((((((9
       Volat( comp%vol+h ) =  sqrt( r1**2 -w**2)
!         if w1 <= 0, it means, h < 0
       if(w .le.  0.) then
          Volat( comp%vol+h ) = - Volat( comp%vol+h ) 
       endif
!))))))))))))
       end
!   ***************************************
      subroutine epbcap(comp, pos, dir, length, icon)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
       integer ir, or, w1, h
       parameter( ir=1, or=2, w1=3, h=4)
!
!          find length to the boundary of 'comp' from 'pos'
!        with direction cos 'dir'
!     'pos' and 'dir' are given in this 'comp' local coordinate
!      in the canonical form.
! 

       type(Component):: comp  ! input. you can extract volume parameters
                               !  by comp.vol(1), etc
       type(epPos)::  pos   ! input.  position.
       type(epDirec)::  dir  ! input. direction cosinse

       real*8  length !  output length cm from pos to the boundary
       integer icon  ! output 0: length obtained. cTrack is inside
                     !        1:  //                        outside
                     !       -1: the line dose not cross the volume
!
       real*8 lout, lin, sp
       integer  jcon1, jcon2

#ifdef DEBUG
       debug = -1
#endif
!          get crossing point with the outer sphere
!(((((((((((
       call kxplsph(pos%x, pos%y, pos%z, dir%x, dir%y, dir%z, 
     *       Volat( comp%vol+or ), lout, jcon1)
!)))))))))
!           jcon =0; pos is inside
!                 1:        outside
!                -1; no cross
       if(jcon1 .ne. -1) then
!           get scaler product for later use
          sp = pos%x*dir%x + pos%y*dir%y +  pos%z*dir%z

!          get crossing point with the inner sphere
!(((((((((((
          call kxplsph(pos%x, pos%y, pos%z, dir%x, dir%y, dir%z,
     *        Volat( comp%vol+ir ), lin, jcon2)
!))))))))))))
          if(jcon1 .eq. 1 .and. jcon2 .eq. -1 ) then
!
!             cross only with outer sphere
!
!(((((((((
             if(lout*dir%z + pos%z .ge. Volat( comp%vol+h ) )  then
!))))))))))
!                    lout should be the solution
                length = lout
                icon = 1
#ifdef DEBUG
                debug = 1
#endif
!((((((((
             elseif( (-lout -2*sp)*dir%z + pos%z  .ge.
     *              Volat( comp%vol+h) ) then
!))))))))
!                should cross at h
!(((((((((
                length =( Volat( comp%vol+h )- pos%z )/dir%z
!))))))))
                icon = 1
#ifdef DEBUG
                debug = 2
#endif
             else
                icon = -1
             endif


          elseif(jcon1 .eq. 1 .and. jcon2 .eq. 1) then

!               penetrate two sphers from  outside
!(((((((((
             if(lout*dir%z + pos%z .ge. Volat( comp%vol+h )) then
!)))))))))
                length = lout
                icon = 1
#ifdef DEBUG
                debug = 3
#endif
             elseif(dir%z .le. 0.) then
                icon = -1
!((((((((((
             elseif((-lout - sp*2)*dir%z + pos%z .lt. 
     *                             Volat( comp%vol+ h)) then
!                      |  this is length to further cross point
                icon =-1
             elseif(lin*dir%z + pos%z .ge. Volat(comp%vol+h)) then
!                cross at h
                length =( Volat(comp%vol+h)- pos%z ) / dir%z
!)))))))))
                icon = 1                
#ifdef DEBUG
                debug = 4
#endif
             elseif((-lin - sp*2)*dir%z + pos%z  .ge.
     &          Volat( comp%vol+h)) then
!                cross with inner sphere
                length = -lin -sp*2
                icon = 1
#ifdef DEBUG
                debug = 5
#endif
             elseif((-lout-2*sp)*dir%z + pos%z .ge. 
     *              Volat( comp%vol+h)) then
!                should cross at h
                length =( Volat( comp%vol+h)- pos%z ) / dir%z
                icon = 1
#ifdef DEBUG
                debug = 6
#endif
             else
                icon = -1
             endif


          elseif(jcon1 .eq. 0 .and. jcon2 .eq. -1) then

!              
!              inside outer sphere and not cross with inner sphere
!
             if(lout*dir%z + pos%z .ge. Volat( comp%vol+h)) then
!                   cross at > h
                if(pos%z .ge. Volat( comp%vol+h)) then
!                     pos > h
                   length = lout
                   icon = 0
#ifdef DEBUG
                   debug = 7
#endif
                else
!                    pos < h; should cross at h                   
                   length =( Volat( comp%vol+h)- pos%z ) / dir%z
                   icon = 1
#ifdef DEBUG
                   debug = 8
#endif
                endif
             else
!                    cross at < h
                if(pos%z .le. Volat( comp%vol+h) ) then
!                       pos < h
                   icon = -1
                else
!                   pos> h; should cross at h
                   length =( Volat( comp%vol+h)- pos%z ) / dir%z
                   icon = 0
#ifdef DEBUG
                   debug = 9
#endif
                endif
             endif
             
          elseif(jcon1 .eq. 0  .and.  jcon2 .eq. 1)  then
! 
!              inside outer sphere and outside inner sphere
!              and crosses the inner spere
!
             if(pos%z .ge. Volat( comp%vol+h)) then
!                     pos > h
                if(lin*dir%z + pos%z .ge. Volat( comp%vol+h)) then
                   length = lin
                   icon = 0
#ifdef DEBUG
                   debug = 10
#endif 
                else
!                    should cross at h
                   length =( Volat( comp%vol+h)- pos%z ) / dir%z
                   icon = 0
#ifdef DEBUG
                   debug = 11
#endif
                endif
             else
!                   pos < h
                if(lout * dir%z + pos%z .lt. Volat( comp%vol+h)) then
                   icon = -1
                elseif( lin*dir%z + pos%z .lt. Volat( comp%vol+h) .and.
     *              ( (-lin -2*sp)*dir%z + pos%z .ge.
     *                     Volat( comp%vol+h)))  then
                   length = (-lin -2*sp)
                   icon = 1
#ifdef DEBUG
                   debug = 12
#endif
                else
!                    should cross at h
                   length =( Volat( comp%vol+h)- pos%z ) / dir%z
                   icon = 1
#ifdef DEBUG
                   debug = 13
#endif
                endif
             endif


          else
!
!                should be inside inner sphere
!
             if(lin*dir%z +  pos%z .ge. Volat( comp%vol+h)) then
!                    cross > h
                length = lin
                icon = 1
#ifdef DEBUG
                debug = 14
#endif
             elseif(lout *dir%z + pos%z  .ge. Volat( comp%vol+h)) then
!                 should cross at h
                length =( Volat( comp%vol+h)- pos%z ) / dir%z
                icon = 1
#ifdef DEBUG
                debug = 15
#endif
             else
                icon = -1
             endif
          endif      
       else
          icon = -1
       endif
       end
!      **********************************
      subroutine epscap(comp, pos, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

       integer ir, or, w1, h
       parameter( ir=1, or=2, w1=3, h=4)

!                judges if pos is inside the comp component      
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord. of ncx-th comp.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside

      real*8 temp
      if(pos%z .lt. Volat( comp%vol+h)) then
         icon = 1
      elseif(pos%z .gt. Volat( comp%vol+or)) then
         icon = 1
      elseif(abs(pos%x) .gt. Volat( comp%vol+or)) then
         icon =1 
      elseif(abs(pos%y) .gt. Volat( comp%vol+or)) then
         icon =1 
      else
         temp = pos%x**2 + pos%y**2 + pos%z**2 
         if(temp .ge. Volat( comp%vol+ir)**2 .and.
     *      temp .le. Volat( comp%vol+or)**2 ) then
            icon = 0
         else
            icon = 1
         endif
      endif
      end
!     **************************************
      subroutine epenvlpcap(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

!
!        give the envloping box of the component
       type(Component)::  comp  ! input component.
       type(epPos)::  org       ! output.  origin of the enveloping box
       type(ep3Vec)::  abc      ! output.  a,b,c of the box

      integer ir, or, w1, h
      parameter( ir=1, or=2, w1=3, h=4)
!
      
      if(Volat( comp%vol+w1) .gt. 0.) then
         org%x = 
     *   - sqrt( Volat( comp%vol+or)**2 - Volat( comp%vol+h)**2)
      else
         org%x = - Volat( comp%vol+or)
      endif
      org%y = org%x
      org%z = Volat( comp%vol+h)
      abc%x = abs(org%x)*2
      abc%y = abc%x
      abc%z = Volat( comp%vol+or) - Volat( comp%vol+h)
      end
!     *************************************
      subroutine epatloccap(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp ! input.
      integer loc(*)   ! output. Volat(comp.vol+loc(i)) is the
                      !  the vlaue of the i-th original
                      !  attribute in the config file.

      integer i

      do i = 1, 4
         loc(i) = i
      enddo

      end

