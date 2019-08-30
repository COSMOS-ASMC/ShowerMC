!     *************************
!     * magage eliptic cylinder
!     ************************
      subroutine eprecyl(comp)
       implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepManager.h"
!
!            read eliptic cylinder configuration
!
       type(Component)::  comp  ! ouptut. here config data is put.

       character*120 msg
 
       integer ra, rb, he
       parameter( ra=1, rb=2, he=3)

       real*8 a, b, h
!
!           read eliptic cylinder configuration data 
!           eliptic cylinder has 3 volume  attributes, 
!           For direction,  its x and y axis. 
!
        call eprpst(comp, 3, 3, 1, 6)

!           check some values
        a = Volat( comp%vol + ra)
        b = Volat( comp%vol + rb)
        h = Volat( comp%vol + he)
        if(a .eq. b) then
           write(msg, *) comp%cn,'-th component: a=', a,
     *    ' = b =', b, ' for ecyl'
           if(MsgLevel-1 .ge. 0) then
              call cerrorMsg(msg, 1)
           endif
        endif
        if(a .le. 0 .or. b .lt. 0. .or.  h .le. 0) then
           write(msg, *) 
     *     comp%cn,'-th component: some attrib <=0 ', a, b, h,
     *    ' for ecyl'
           call cerrorMsg(msg, 0)
        endif
      end
!   ***************************************
      subroutine epbecyl(comp, pos, dir, length, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"
#include "ZepDirec.h"
#include "Zepdebug.h"
!
 
       integer ra, rb, he
       parameter( ra=1, rb=2, he=3)

       real*8 a, b, h
!    
!
!        find length to the boundary from the pos with dir
!        
       type(Component):: comp  ! input.  boundary to this comp. is examined
       type(epPos)::  pos    ! input. local coordinate position
       type(epDirec)::  dir   ! input. local coord. direction cos.
      real*8  length  ! output length cm
      integer icon    ! output
                      ! 0: length obtained. cTrack is inside
                      ! 1:  //                        outside
                      !-1: the line dose not cross the volume


      real*8  d, inout, alpha, beta,  crossh, sqrtd, crossh2


#ifdef  DEBUG
      debug=-1
#endif
      a = Volat( comp%vol + ra)
      b = Volat( comp%vol + rb)
      h = Volat( comp%vol + he)
      alpha =  (dir%x/a)**2 + (dir%y/b)**2
      if(alpha .eq. 0.) then
!          in some case, |dir.z| may not be exactly 1
!          so correct it
         if(dir%z .lt. 0.) then
            dir%z=-1.d0
         else
            dir%z= 1.0
         endif
      endif
      beta = pos%x*dir%x/a/a + pos%y*dir%y/b/b
      inout =(pos%x/a)**2 + (pos%y/b)**2 - 1.d0
      d = beta**2 - alpha * inout
      if(d .lt. 0.d0 .or. alpha .eq. 0.d0) then
         if(abs(dir%z) .eq. 1.d0 .and. inout .le. 0.d0) then
            if(pos%z .ge. h) then
               if( dir%z .lt. 0. ) then
                  length = pos%z - h
                  icon = 1
#ifdef DEBUG                
                  debug = 8
#endif
               else
                  icon = -1
               endif
            elseif(pos%z .lt. h .and. pos%z .gt. 0. ) then
               if(dir%z .gt .0.) then
                  length = h - pos%z
                  icon = 0
#ifdef DEBUG                
                  debug = 9
#endif
               else
                  length = pos%z
                  icon = 0
#ifdef DEBUG
                  debug = 10
#endif
               endif
            else
!                    pos.z <=0
               if(dir%z .gt. 0.) then
                  length = -pos%z
                  icon = 1
#ifdef DEBUG
                  debug = 11
#endif
               else
                  icon = -1
               endif
            endif
         else
            icon = -1
         endif
      elseif(inout .le. 0.d0) then
!               posz is in ellipse
         length = (-beta + sqrt(d))/alpha
         if(pos%z .ge. 0.d0 .and.  pos%z .le. h) then
            crossh = pos%z +  length*dir%z
            if( crossh .gt. h ) then
               icon = 0
               length = (h-pos%z)/dir%z
#ifdef  DEBUG
               debug=1
#endif
            elseif(crossh .lt. 0.) then
               length = -pos%z/dir%z
               icon = 0
#ifdef  DEBUG
               debug=7
#endif
            else
               icon = 0
#ifdef  DEBUG
               debug=0
#endif
            endif
         elseif(pos%z .gt. h) then
            if(pos%z + length*dir%z .le. h) then
!                should cross at h                
               length = (h - pos%z)/ dir%z
               icon  = 1
#ifdef  DEBUG
               debug=2
#endif
            else
               icon = -1
            endif
         else
!            pos.z < 0
            if(pos%z + length*dir%z .ge. 0.d0) then
               length = -pos%z/dir%z
               icon = 1
#ifdef  DEBUG
               debug=3
#endif
            else
               icon =-1
            endif
         endif
      else
!          point is outside cyl
         sqrtd = sqrt(d)
         length = (-beta - sqrtd)/alpha
         if( length .lt. 0.) then
            icon =-1
            return     !***********
         endif
         crossh = pos%z + length *dir%z 
         if(crossh .ge. 0. .and. crossh .le. h) then
            icon = 1
#ifdef  DEBUG
            debug=4
#endif
         else
!                 further point
            length = (-beta + sqrtd)/alpha 
            if(length .lt. 0.) then
               icon =-1  
               return  ! ***********
            endif
            crossh2 = pos%z + length*dir%z
            if(crossh .gt. h .and. crossh2 .lt. h) then
               length = (h-pos%z)/dir%z
               icon = 1
#ifdef  DEBUG
               debug = 5
#endif
            elseif(crossh .lt. 0. .and. crossh2 .gt. 0.) then
               length = -pos%z/dir%z
               icon = 1
#ifdef  DEBUG
               debug = 6
#endif
            else
               icon = -1 
            endif
         endif
      endif
      end
!      **********************************
      subroutine epsecyl(comp, pos, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


!                judges if pos is inside the comp component      
       type(Component)::  comp !input component
       type(epPos)::  pos  ! input. position in  local coord. of ncx-th comp.
      integer icon  ! output. 0--> pos is inside
                    !         1-->        outside
 
       integer ra, rb, he
       parameter( ra=1, rb=2, he=3)




      if(pos%z .gt. Volat( comp%vol + he)) then
         icon = 1
      elseif(pos%z .lt. 0.) then
         icon = 1
      elseif( (pos%x/Volat( comp%vol + ra))**2 +
     *      (pos%y/Volat( comp%vol + rb))**2 .gt. 1.d0) then
         icon =1 
      else
         icon =0
      endif

      end
!     **************************************
      subroutine epenvlpecyl(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

!
!        give the envloping box of the component
       type(Component)::  comp  ! input component.
       type(epPos)::  org       ! output.  origin of the enveloping box
       type(ep3Vec)::  abc      ! output.  a,b,c of the box
 
       integer ra, rb, he
       parameter( ra=1, rb=2, he=3)


      org%x =  -Volat( comp%vol + ra)
      org%y =  -Volat( comp%vol + rb) 
      org%z =  0.
      abc%x =  Volat( comp%vol + ra)*2
      abc%y =  Volat( comp%vol + rb)*2
      abc%z =  Volat( comp%vol + he) 
      end
         
!     **************************************
      subroutine epatlocecyl(comp, loc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"

       type(Component)::  comp  ! input component.
      integer loc(3)      

       loc(1) = 1
       loc(2) = 2
       loc(3) = 3
      end
         


