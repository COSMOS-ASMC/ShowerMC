      subroutine ctransCoord2(sys, a, b)
!        sys: character string. input target coordinate system.
!          a: /coord/  input.
!          b: /coord/  outupt.  b containes transformed coordinate.
!                               b can be a.
! a.sys  \sys        'xyz'    'llh'   'sph'   
!
!  'xyz'               o        o       o     
!  'llh'               o        o       o   
!  'sph'               o        o       o
! 
!
      implicit none

#include  "Zcoord.h"


      type(coord)::a
      type(coord)::b
      character*(*) sys
      character*70 msg
!
      if(a%sys .eq. sys) then
!             already in the objective system
          b = a
      else
         if(a%sys .eq. 'xyz') then     
            if(sys .eq. 'llh') then  
               call ceCent2llh(a, b)
            elseif(sys .eq. 'sph') then 
               call ceCent2sph(a, b)
            else   
               write(msg, *) 'error sys=',sys,' transCoord2'
               call cerrorMsg(msg, 0)
            endif
         elseif(a%sys .eq. 'llh') then
            if(sys .eq.  'xyz') then
                call cllh2eCent(a, b)
            elseif(sys .eq. 'sph') then  
               call cllh2sph(a, b)
            else
               write(msg, *) 'error sys=',sys,' transCoord2'
               call cerrorMsg(msg, 0)
            endif
         elseif(a%sys .eq. 'sph') then   ! from 'sph'
            if(sys  .eq.  'xyz') then    
               call csph2eCent(a, b)                ! to 'xyz'
            elseif(sys .eq. 'llh') then
               call csph2llh(a, b)                 ! to 'llh'
            else   
               write(msg, *) 'error a%sys=',a%sys,' transCoord2'
               call cerrorMsg(msg, 0) 
            endif
         else
               write(msg, *) 'error a%sys=',a%sys,' transCoord2'
               call cerrorMsg(msg, 0)
         endif  
       endif
       end
