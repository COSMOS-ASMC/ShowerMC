      subroutine epbndry2(cnx, el, icon)
      use moddebug
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
!
!          get the length to the nearest crossing point of cTrack
!          and the specified component
!
        integer cnx   ! input. cnx-th component is examined
        real*8  el    ! output. length to the boundary of the cnx-th  comp
                      !         from cTrack.
        integer icon  ! output. 0, x,y,z is inside
                      !         1  outside,
                      !        -1  no Xssing

       type(epPos)::  posw, posl
       type(epDirec)::  dirw, dirl


        character*80 msg
        character(len=MAX_STRUCCHR):: epparaphrase, tempph
        character(len=MAX_STRUCCHR):: basename
        integer uscl
         if(Cn .eq. cnx) then
!            no need to conv. to local.coor.
           posl =cTrack%pos
           dirl =cTrack%w 
        else
           call epl2w(Cn, cTrack%pos, posw)
           call epw2l(cnx, posw, posl)
           call epl2wd(Cn, cTrack%w, dirw)
           call epw2ld(cnx, dirw, dirl)
        endif

        call epGetBaseStrucName(Det%cmp(cnx)%struc, basename)
        
        if(basename .eq. 'box') then
           call epbbox(Det%cmp(cnx), posl, dirl, el, icon)
        elseif(basename .eq. 'cyl') then
           call epbcyl(Det%cmp(cnx), posl, dirl, el, icon)
        elseif(basename .eq. 'pipe') then
           call epbpip(Det%cmp(cnx), posl, dirl, el, icon)
        elseif(basename .eq. 'prism' ) then
           call epbprs(Det%cmp(cnx), posl, dirl,  el, icon)
        elseif(basename .eq. 'sphere') then
           call epbsph(Det%cmp(cnx), posl, dirl, el, icon)
        elseif(Det%cmp(cnx)%struc(1:4) .eq. 'new-') then
           call epbNew(Det%cmp(cnx), posl, dirl, el, icon)
        else
           call epseeUnderScore(Det%cmp(cnx)%struc, uscl)
           tempph = epparaphrase(Det%cmp(cnx)%struc(1:uscl))
           if(tempph(1:4) .eq. 'new-') then
              call epbNew(Det%cmp(cnx), posl, dirl, el, icon)
           else
!                   
              write(msg, *) ' strange structure=',
     *             Det%cmp(cnx)%struc
              call cerrorMsg(msg, 0)
           endif
        endif
        end
