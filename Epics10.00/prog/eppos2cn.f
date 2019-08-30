       subroutine eppos2cn(ncx, aTrack, nc)
       implicit none

!          find  component # of given track which is
!          in or near ncx-th component. 
!          if ncx = 0, 
!               aTrack.pos is assumed to be in world coord.
!          else
!               aTrack.pos is assumed to be in local coord.
!               of the ncx-th comp.
!
#include  "ZepTrack.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"

       integer ncx
       type(epTrack)::  aTrack  ! input. ask  where is this track
       integer nc               ! output.  found comp. #, if
                        ! the point is outside of any comp.
                        ! nc=Det.nct+1

       character*70 msg
!
!!       if(form .eq. 'box') then
!!c                   only boxes with the same a, b 
!!          call epsbox(ncx, aTrack, nc)
!!       elseif(form .eq. 'cyl') then
!!          call epscyl(ncx, aTrack, nc)
!!       elseif(form .eq. 'pipe') then
!!          call epspip(ncx, aTrack, nc)
!!       elseif(form .eq.  'mix') then
          call epsall(ncx, aTrack%pos,  nc)
!       else
!          write(msg,*)
!     *         ' eppos2cn for form=',form, ' not supported'
!          call cerrorMsg(msg, 0)
!       endif
      end
      subroutine epsInside(pos, nci,  nco)
      implicit none
#include  "ZepPos.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       type(epPos)::  pos  !  input. see if pos is in some matreska.
!                                   pos is in W.C
      integer   nci  !  input.  component number  which  includes pos,
                     !  and is Matreska
      integer   nco  !  output.  Matreska's  component number which
                     !          includes pos.  if 0, there is no
                     !          child matreska which includes pos.

       type(epPos)::  posl

      integer i, icon
      do i  = 1, Det%cmp(nci)%NMatreska
!(((((((((((((((
!         call epw2l(Det.cmp(nci).Contains(i), pos, posl)
!         call epsOne(Det.cmp(nci).Contains(i), posl, icon)
         call epw2l(CnArea( Det%cmp(nci)%Contains+i), pos, posl)
         call epsOne(CnArea( Det%cmp(nci)%Contains+i), posl, icon)
!)))))))))))))
         if(icon .eq. 0) then
!((((((((((((((
!            nco  = Det.cmp(nci).Contains(i)
            nco  = CnArea( Det%cmp(nci)%Contains+i )
!))))))))))
            goto 10
         endif
      enddo
      nco =  0
 10   continue
      end


      subroutine epsOne(ncx, pos, icon)
      implicit none
#include  "ZepPos.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"

!          see if pos is in ncx-th comp.
!          pos is in l.c of  ncx
!          if  yes, icon = 0 else  icon =1
!        This simply judge pos is inside  of ncx-th comp.
!       Even if yes, the point may be in an inner componennt.
!       or not included in a union part.
!       So that it must be examined further.
!
       type(epPos)::  pos
      
      character*80 msg
      character(len=MAX_STRUCCHR)::  epparaphrase, tempph

      integer  ncx,  i, icon
      integer uscl
!     
      i = ncx

      if(Det%cmp(i)%struc  .eq. 'box' .or.
     *   Det%cmp(i)%struc(1:4)  .eq. 'box_' ) then
         call epssbox(i, pos, icon)

      elseif(Det%cmp(i)%struc .eq. 'cyl' .or.
     *       Det%cmp(i)%struc(1:4) .eq. 'cyl_') then
         call epscyl1(i, pos, icon)

      elseif(Det%cmp(i)%struc .eq. 'pipe' .or.
     *       Det%cmp(i)%struc(1:5) .eq. 'pipe_' ) then
         call epspip1(i, pos, icon)

      elseif(Det%cmp(i)%struc .eq. 'sphere' .or.
     *       Det%cmp(i)%struc(1:7) .eq. 'sphere_' ) then
         call epsssph(i, pos, icon)

      elseif(Det%cmp(i)%struc .eq. 'prism' .or.
     *      Det%cmp(i)%struc(1:6) .eq. 'prism_' ) then   
         call epsprism(i, pos, icon)
      elseif(Det%cmp(i)%struc(1:4) .eq. 'new-') then
         call epsNew(Det%cmp(i), pos, icon)
      else
         call epseeUnderScore(Det%cmp(i)%struc, uscl)
         tempph = epparaphrase(Det%cmp(i)%struc(1:uscl))
         if(tempph(1:4) .eq. 'new-') then
            call epsNew(Det%cmp(i), pos, icon)
         else
            write(msg, *) ' struc=', trim(Det%cmp(i)%struc),
     *           '  for ',i,'-th comp. invalid'
            call cerrorMsg(msg, 0)
         endif
      endif
      end
      subroutine epsPCont(ncx, pos, icon)
      implicit none
#include  "ZepPos.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!          see if pos, in ncx-th comp.,
!          is included in a parent which partially
!          includes ncx-th comp.  pos is in ncx-local.
!          if  yes, icon = 0 else  icon =1
!          If there is no parent, icon =0
!
       type(epPos)::  pos

      integer  ncx,  icon

      character*80 msg
      character(len=MAX_STRUCCHR)::  epparaphrase, tempph

      integer i, j, ii

       type(epPos)::  posw, posl
!           see if parent matreska  which partially
!           contains this ncx-th comp.
!           If the point is not in any of such a parent
!           matreska, icon = 1 will result
!           so the point is considered not included in
!           the current component.
      character(len=MAX_STRUCCHR):: basename

      i = ncx
      if(Det%cmp(i)%NPContainer .gt. 0) then
         call epl2w(ncx, pos, posw)
         do j =1, Det%cmp(i)%NPContainer
!                see if pos is inside the  parent matreshka which is ii-th
!                comp.
!((((((((((((((
!            ii = Det.cmp(i).PContained(j)
            ii = CnArea( Det%cmp(i)%PContained+j )
!))))))))))
            call epw2l(ii, posw,  posl)
            call epGetBaseStrucName(Det%cmp(ii)%struc, basename)

            if(basename == 'box') then
               call epssbox(ii, posl, icon)
            elseif(basename == 'cyl') then
               call epscyl1(ii, posl, icon)
            elseif(basename == 'pipe') then
               call epspip1(ii, posl, icon)
            elseif( basename == 'sphere') then
               call epsssph(ii, posl, icon)
            elseif( basename == 'prism' ) then
               call epsprism(ii, posl,  icon)
            elseif(Det%cmp(ii)%struc(1:4) .eq. 'new-') then
               call epsNew(Det%cmp(ii),  posl, icon)
            else
               tempph = epparaphrase(Det%cmp(ii)%struc)
               if(tempph(1:4) .eq. 'new-') then
                  call epsNew(Det%cmp(ii),  posl, icon)
               else
                  write(msg, *)
     *                ' struc=', Det%cmp(i)%struc,'  for ',i,
     *                '-th comp. invalid'
                  call cerrorMsg(msg, 0)
               endif
            endif
            if(icon .eq. 0) goto 100
         enddo
      else
         icon = 0
      endif
 100  continue
      end
!             
      subroutine  epsMatreska(ncx, pos, nc)
      implicit none
#include  "ZepPos.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!          see if pos  in ncx-th comp.
!          is included in one of the child Matreskas.
!          pos is in ncx-local.
!          if some matreska is found, nc will get
!          that comp. # or nc = ncx.
!       
!      
       type(epPos)::  pos
      integer  ncx,    nc
      integer  nco
       type(epPos)::  posw

      nc = ncx
      if(Det%cmp(nc)%NMatreska .gt. 0) then

         call epl2w(nc, pos, posw)

         do  while(Det%cmp(nc)%NMatreska .gt. 0)
!           find inner most matreska
!           nc-th comp.is matreska,  search
!           if posw, is inside  matreska
!              posw in w.c
            call epsInside(posw, nc, nco)

!                       nco =0 means there is no matreska to contain
!                       pos
            if(nco .eq. 0) goto 300
!                 nco-th comp. contains pos, see if more matreska
            nc = nco
         enddo
      endif
 300  continue
      end

      subroutine epsall(ncx, pos, nc)
!          ncx=0--> pos is in w.c
!         else      in l.c of ncx.
      implicit none
#include  "ZepPos.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"

       type(epPos)::  pos
      integer  ncx, nc
      integer i, icon, jcon
       type(epPos)::  posw, posl

!        assume  we may search only box or sphere_w
!        which is uncontained first.
      if(ncx .eq. 0)  then
         posw = pos
      else
!             search firstly  ncx-comp for which no conversion is
!             needed.
!             pos is in ncx-th comp.
         call epsOne(ncx, pos, icon)

         if(icon .eq. 0) then
!              yes, so examine matreska
            call epsMatreska(ncx, pos, nc)
!/////////////////
!            call debugpos('epspos2cn matreska', ncx, pos, nc)
!//////////////////
!              if not matreska, nc=ncx. if pos in inner comp.
!              nc becomes the inner most one.
            if(ncx .ne. nc) then
               call epl2w(ncx, pos, posw)
               call epw2l(nc,  posw, posl)
            else
               posl = pos
            endif
!            call epsPCont(nc, pos, jcon)
            call epsPCont(nc, posl, jcon)
! /////////////
!            call debugpos('epspos2cn epsPCont',nc,pos,jcon)
!/////////////////
            if(jcon .ne. 0) then
!               nc is  partially contained by other comp.
!               but  pos is not in that part. i.e., outside
               call epl2w(ncx, pos, posw)
               goto 5
            endif
            goto 10
         else
!                prepare to search in other comp. covert
!                pos into posw
            call epl2w(ncx, pos, posw)
         endif
      endif
 5    continue
!                 -----------------
!       search other comp.  for the fast search, 
!       first examine uncontained component  only.
!!!!!!! v9.14 except for small detectors consisting of same
!             area, only the world is uncontained, so
!             we may search from the last one
!      do i = ncx+1,  Det.nct
      do i =  Det%nct, 1, -1
!         if( (Det.cmp(i).struc(1:3) .eq. 'box'  .or.
!     *        Det.cmp(i).struc .eq. 'sphere_w')  
!     *     .and.
!     *        Det.cmp(i).NContainer .eq. 0) then
         !!!!!!!!v9.14 only world is with NContainer= 0
         !!!! so skip non-world comp. 
!//////////
!         write(0,*) ' i=',i,  ' Det.cmp(i).NContainer=', 
!     *        Det.cmp(i).NContainer
!         write(0,*) ' i=',i,  ' Det.cmp(i).NPContainer=', 
!     *        Det.cmp(i).NPContainer
!         write(0,*) ' i=',1,  ' Det.cmp(1).NPContainer=', 
!     *        Det.cmp(1).NPContainer
!////////////
         if( i < Det%nct .and. Det%nworld == 1 ) then
            nc= Det%nct+1
!////////////
!            write(0,*) ' i=',i,' < Det.nct && nworld=1'
!            write(0,*) 'so goto 10'
!//////////////
            goto 10
         endif
         !!!!!!!!
         if(Det%cmp(i)%NContainer .eq. 0 .and.   !! added
     *      Det%cmp(i)%NPContainer == 0 ) then
            call epw2l(i, posw, posl)
!//////////
!            write(0,*) ' posw=',posw.x,posw.y, posw.z
!//////////////////
!            write(0,*) ' converted to posl for i=',i
!            write(0,*) 'posl=',posl.x,posl.y, posl.z
!            write(0,*) '  posl is examined by epsOne'
!/////////////////
            call epsOne(i, posl, icon)
!///////////////
!            call debugpos('eppos2cn;epsOne', i, posl,icon)
!//////////////////
            if(icon .eq.  0) then
               call epsMatreska(i, posl, nc)
!///////////////
!             call debugpos('eppos2cn;epsOne;matre', i, posl,icon)
!//////////////////
               if(nc .ne. i) then
                  call epw2l(nc, posw, posl)
               endif
               call epsPCont(nc, posl, jcon)
!///////////////
!            call debugpos('eppos2cn;epsOne;PCont', nc, posl,jcon)
!//////////////////
               if(jcon .ne. 0) then
                  nc = Det%nct + 1
               endif
               goto 10
            endif
         endif
      enddo
!      do i = ncx -1, 1, -1
!c         if( (Det.cmp(i).struc(1:3) .eq. 'box'  .or.
!c     *        Det.cmp(i).struc .eq. 'sphere_w')  
!c     *     .and.
!c     *        Det.cmp(i).NContainer .eq. 0) then
!         if(Det.cmp(i).NContainer .eq. 0) then
!
!            call epw2l(i, posw, posl)
!            call epsOne(i, posl, icon)
!c///////////////
!c            call debugpos('eppos2cn;epsOne;inv', i, posl,icon)
!c//////////////////
!
!            if(icon .eq. 0) then
!               call epsMatreska(i, posl, nc)
!c///////////////
!c            call debugpos('eppos2cn;matre;inv', i, posl,nc)
!c//////////////////
!
!               if(i .ne. nc) then
!                  call epw2l(nc, posw, posl)
!               endif
!               call epsPCont(nc, posl, jcon)
!c///////////////
!c            call debugpos('eppos2cn;mPCont;inv', nc, posl, jcon) 
!c//////////////////
!
!               if(jcon .ne. 0) then
!                  nc = Det.nct + 1
!               endif
!               goto 10
!            endif
!         endif
!      enddo
      nc= Det%nct+1

 10   continue
!///////////
!      write(0,*) ' exiting epsall: nc=',nc
!////////////
      end
!     ***************************************
      subroutine epsbox(ncx, aTrack, nc)
      implicit none
!          find  component # of position pos which is
!          in or near ncx-th component. 
!          if ncx =0, xx is world coord.
!          else ncx-local coord.
#include  "ZepTrackv.h"
#include  "Zcnfig.h"

       type(epTrack)::  aTrack
       type(epPos)::  pos
      integer ncx, nc

!
      pos = aTrack%pos
      if(ncx .le. 0) then
         call epsbxa(pos, nc)
      elseif(ncx .le. Det%nct) then
         if(0.d0 .le. pos%z      .and.
     *        pos%z .le. Volat( Det%cmp(ncx)%vol + boxc)) then
            if(0.d0 .le. pos%x .and.
     *           pos%x .le. Volat( Det%cmp(ncx)%vol + boxa) ) then
               if(0.d0 .le. pos%y .and.
     *              pos%y .le. Volat( Det%cmp(ncx)%vol + boxb)) then
                  nc=ncx
               else
                  nc=Det%nct+1
               endif
            else
               nc=Det%nct+1
            endif
         else
            if(aTrack%w%z .gt. 0.) then
               call epsbxu(ncx, pos, nc)
            else
               call epsbxd(ncx, pos, nc)
            endif
         endif
      else
         nc=Det%nct+1
      endif
      end
!!c      *****************
!!      subroutine epscyl(ncx, aTrack, nc)
!!      implicit none
!!#include  "ZepTrackv.h"
!!#include  "Zcnfig.h"
!!
!!      record /epTrack/ aTrack
!!      record /epPos/ pos
!!      integer ncx, nc
!!      record /epPos/ post
!!
!!c            ncx=0 --> pos is w.c
!!c               else ncx-local
!!      post = aTrack.pos
!!
!!      if(ncx .le. 0) then
!!         call epscya(post, nc)
!!      elseif(ncx .le. Det.nct) then
!!C&&&&&&&
!!         call epv2c_cyl(Det.cmp(ncx), post, pos)
!!c&&&&&
!!         if(0.d0 .le. pos.z      .and.
!!     *        pos.z .le. Volat( Det.cmp(ncx).vol + cylh)) then
!!            if(pos.x**2 + pos.y**2 .le.
!!     *           Volat( Det.cmp(ncx).vol + cylr)**2) then
!!               nc=ncx
!!            else
!!               nc = Det.nct+1
!!            endif
!!         else
!!            if(aTrack.w.z .gt. 0.) then
!!               call epscyu(ncx, pos, nc)
!!            else
!!               call epscyd(ncx, pos, nc)
!!            endif
!!         endif
!!      else
!!         nc = Det.nct+1
!!      endif
!!      end
! ******************************* search only ncx-th cyl
      subroutine epscyl1(ncx, posi, icon)
        implicit none
#include "ZepTrackv.h"
#include "Zcnfig.h"

       type(epPos)::  posi
        integer ncx, icon
!    &&&&&&&&&&
       type(epPos)::  pos


!          pos in local-c.
!  &&&&         convert to canonial one
        call epv2c_cyl(Det%cmp(ncx), posi, pos)
!  &&&&&&
        if(0.d0 .le. pos%z      .and.
     *       pos%z .le. Volat( Det%cmp(ncx)%vol + cylh)) then
           if(pos%x**2 + pos%y**2 .le.
     *          Volat( Det%cmp(ncx)%vol + cylr)**2) then
              icon  =0
           else
              icon = 1
           endif
        else
           icon = 1
        endif
        end
!!c ****************************
!!       subroutine epspip(ncx, aTrack, nc)
!!       implicit none
!!#include "ZepTrackv.h"
!!#include "Zcnfig.h"
!!
!!       record /epTrack/ aTrack
!!       integer ncx,  nc
!!        
!!       real*8 temp
!!       record /epPos/ pos
!!c          ncx =0  --> xx in w.c
!!c           else     in l.c
!!
!!       pos = aTrack.pos
!!
!!       if(ncx .le. 0) then
!!          call epspipa(pos, nc)
!!       elseif(ncx .le. Det.nct) then
!!          if(0.d0 .le. pos.z      .and.
!!     *         pos.z .le. Volat( Det.cmp(ncx).vol + pipeh)) then
!!             temp = pos.x**2 + pos.y**2 
!!             if( temp .le. Volat( Det.cmp(ncx).vol + pipeor)**2 .and.
!!     *          temp .ge. Volat( Det.cmp(ncx).vol + pipeir)**2 ) then
!!                nc=ncx
!!             else
!!                nc = Det.nct+1
!!             endif
!!          else
!!             if(aTrack.w.z .gt. 0.) then
!!                call epspipu(ncx, pos, nc)
!!             else
!!                call epspipd(ncx, pos, nc)
!!             endif
!!          endif
!!       else
!!          nc = Det.nct+1
!!       endif
!!       end
       subroutine epsbxu(ncx, pos, nc)
       implicit none
!              search upward  when form=box
!       xx, yy, zz is assumed to be in ncx-local
!
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
          
       integer ncx, nc
       type(epPos)::  pos
       
       integer i
       type(epPos)::  posw      ! world coord.
       type(epPos)::  posl      !  local. coord.

       call epl2w(ncx, pos, posw)
       do   i=max(ncx-1,1), Det%nct
          if(i .eq. ncx) goto 100
          call epw2l(i, posw, posl)
          if(0.d0 .le. posl%z  .and.
     *         posl%z .le.  Volat( Det%cmp(i)%vol + boxc)) then
             if(0.d0 .le. posl%x .and.
     *            posl%x .le.  Volat( Det%cmp(i)%vol + boxa)) then
                if(0.d0 .le. posl%y .and.
     *               posl%y .le. Volat( Det%cmp(i)%vol + boxb)) then
                   nc=i
                   goto 101
                else
                   nc = Det%nct+1
                   goto 101
                endif
             else
                nc = Det%nct+1
                goto 101
             endif
          elseif(0.d0 .gt. posl%z) then
             nc = Det%nct+1
             goto 101
          endif
 100      continue
       enddo
       nc= Det%nct+1
 101   continue
       end 
!
!       ************
        subroutine epsbxd(ncx, pos, nc)
!       ************
       implicit none
!              search upward  when form=box
!       xx, yy, zz is assumed to be in ncx-local
!
#include  "ZepPos.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
          
           integer ncx, nc
       type(epPos)::  pos

           integer i
       type(epPos)::  posw ! world coord.
       type(epPos)::  posl !  local. coord.

!             search downword
        call epl2w(ncx, pos, posw)
        do   i=min(ncx+1,Det%nct), 1, -1
           if(i .eq. ncx) goto 200
           call epw2l(i, posw, posl)
           if(0.d0 .le. posl%z  .and.
     *          posl%z .le.  Volat( Det%cmp(i)%vol + boxc)) then
              if(0.d0 .le. posl%x .and.
     *             posl%x .le. Volat( Det%cmp(i)%vol + boxa))  then
                 if(0.d0 .le. posl%y .and.
     *                posl%y .le. Volat( Det%cmp(i)%vol + boxb))then
                    nc=i
                    goto 201
                 else
                    nc = Det%nct+1
                    goto 201
                 endif
              else
                 nc = Det%nct+1
                 goto 201
              endif
           elseif(0.d0 .lt. posl%z) then
              nc = Det%nct+1
              goto 201
           endif
 200       continue 
        enddo
        nc = Det%nct+1
 201    continue
       end
!!       subroutine epscyu(ncx, pos, nc)
!!c            pos in l.c
!!       implicit none
!!#include  "ZepTrackv.h"
!!#include  "Zcnfig.h"
!!       record /epPos/ pos
!!       integer ncx, nc
!!       record /epPos/ posw, posl
!!       integer i
!!
!!       call epl2w(ncx, pos, posw)
!!       do   i=ncx+1, Det.nct
!!          call epw2l(i, posw, posl)
!!          if(0.d0 .le. posl.z  .and.
!!     *         posl.z .le.  Volat( Det.cmp(i).vol + cylh)) then
!!             if(posl.x**2 + posl.y**2 
!!     *            .le. Volat( Det.cmp(i).vol + cylr)**2)   then
!!                nc=i
!!                goto 101
!!             else
!!                nc = Det.nct+1
!!                goto 101
!!             endif
!!          endif
!!       enddo
!!       nc = Det.nct+1
!! 101   continue
!!       end
!
!!c       ************
!!      subroutine epscyd(ncx, pos, nc)
!!c       ************
!!c            pos in l.c
!!       implicit none
!!#include  "ZepTrackv.h"
!!#include  "Zcnfig.h"
!!       record /epPos/ pos
!!       integer ncx, nc
!!       record /epPos/ posw, posl
!!       integer i
!!
!!c             search downward
!!        call epl2w(ncx, pos, posw)
!!        do   i=ncx-1, 1, -1
!!           call epw2l(i, posw, posl)
!!           if(0.d0 .le. posl.z  .and.
!!     *              posl.z .le.  Volat( Det.cmp(i).vol + cylh)) then
!!              if(posl.x**2 + posl.y**2 
!!     *          .le. Volat( Det.cmp(i).vol + cylr)**2)  then
!!                 nc=i
!!                 goto 201
!!              else
!!                 nc = Det.nct+1
!!                 goto 201
!!              endif
!!           endif
!!        enddo
!!        nc = Det.nct+1
!! 201    continue
!!       end
!!!      ***********************************
!!       subroutine epspipu(ncx, pos, nc)
!!c         pos in l.c
!!       implicit none
!!#include  "ZepTrackv.h"
!!#include  "Zcnfig.h"
!!       record /epPos/pos
!!       integer ncx, nc
!!       real*8 temp
!!       integer i
!!       record /epPos/ posw, posl
!!       
!!       call epl2w(ncx, pos, posw)
!!       do   i=ncx+1, Det.nct
!!          call epw2l(i, posw, posl)
!!          if(0.d0 .le. posl.z  .and.
!!     *        posl.z .le.  Volat( Det.cmp(i).vol + pipeh)) then
!!             temp = posl.x**2 +  posl.y**2
!!             if(temp .le. Volat( Det.cmp(i).vol + pipeor)**2 .and.
!!     *            temp .ge. Volat( Det.cmp(i).vol + pipeir)**2  )   then
!!                nc=i
!!                goto 101
!!             else
!!                nc = Det.nct+1
!!                goto 101
!!             endif
!!          endif
!!       enddo
!!       nc = Det.nct+1
!! 101   continue
!!       end
!!c
!!c      ************
!!       subroutine  epspipd(ncx, pos, nc)
!!c      ************
!!       implicit none
!!#include  "ZepTrackv.h"
!!#include  "Zcnfig.h"
!!       record /epPos/pos
!!       integer ncx, nc
!!       real*8 temp
!!       integer i
!!       record /epPos/ posw, posl
!!
!!
!!c             search downward
!!        call epl2w(ncx, pos, posw)
!!        do   i=ncx-1, 1, -1
!!           call epw2l(i, posw, posl)
!!           if(0.d0 .le. posl.z  .and.
!!     *          posl.z .le. Volat( Det.cmp(i).vol + pipeh)) then
!!              temp = posl.x**2 + posl.y**2 
!!              if( temp .le. Volat( Det.cmp(i).vol + pipeor)**2 .and.
!!     *            temp .ge. Volat( Det.cmp(i).vol + pipeir)**2)  then
!!                 nc=i
!!                 goto 201
!!              else
!!                 nc =  Det.nct+1
!!                 goto 201
!!              endif
!!           endif
!!        enddo
!!        nc = Det.nct+1
!! 201    continue
!!       end
!      ********************************
       subroutine  epspipa(pos, nc)
!      ********************************
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos):: pos
       integer nc
       real*8 temp
       integer i
       type(epPos)::  posl, post

!             search all
       do   i=1, Det%nct
          call epw2l(i, pos, post)
!&&&&&
          call epv2c_pipe(Det%cmp(i), post, posl)
!&&&&&

          if(0.d0 .le. posl%z  .and. 
     *        posl%z .le. Volat( Det%cmp(i)%vol + pipeh)) then
             temp = posl%x**2 + posl%y**2 
             if(temp .le. Volat( Det%cmp(i)%vol + pipeor)**2 .and.
     *          temp .ge. Volat( Det%cmp(i)%vol + pipeir)**2) then       
                nc=i
                goto 301
             endif
          endif
       enddo
       nc =  Det%nct+1
 301   continue
       end
       subroutine epsbxa(pos, nc)
       implicit none
!              search all component( the same a, b assumed)
!              pos is in the world coord.
!       ************
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos
       integer nc
       
       type(epPos)::  posl        
       integer i


!       search all
       do   i=1, Det%nct
          call epw2l(i, pos, posl)
          if(0.d0 .le. posl%z  .and.
     *         Volat( Det%cmp(i)%vol + boxc) .ge. posl%z) then
             if(0.d0 .le. posl%x .and.
     *            Volat( Det%cmp(i)%vol + boxa) .ge. posl%x) then
                if(0.d0 .le. posl%y .and.
     *               Volat( Det%cmp(i)%vol + boxb) .ge. posl%y) then
                   nc=i
                   goto 301
                else
                   nc = Det%nct+1
                   goto 301
                endif
             else
                nc = Det%nct+1
                goto 301
             endif
          endif
       enddo
       nc = Det%nct+1
 301   continue
       end
       subroutine epscya(pos, nc)
!             pos is in w.c
       implicit none
!       ************
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos
       integer nc
       type(epPos)::  posl, post
       integer i

!             search all
       do   i=1, Det%nct
          call epw2l(i, pos, post)
!&&&&&
          call epv2c_cyl(Det%cmp(i), post, posl)
!&&&&&
          if(0.d0 .le. posl%z  .and. 
     *         Volat( Det%cmp(i)%vol + cylh) .ge. posl%z) then
             if( posl%x**2 + posl%y**2 .le.
     *            Volat( Det%cmp(i)%vol + cylr)**2)  then
                nc=i
                goto 301
             else
                nc =  Det%nct+1
                goto 301
             endif
          endif
       enddo
       nc = Det%nct+1
 301   continue
       end
       subroutine epspip1(ncx, posi, icon)
!         search only ncx-th pipe
!         xx .. is assumed to be in local coord.
       implicit none
!       ************
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  posi
       integer  icon, ncx
       real*8 temp

       type(epPos)::  pos

!  &&&&         convert to canonial one
        call epv2c_pipe(Det%cmp(ncx), posi, pos)
!  &&&&&&


       if(0.d0 .le. pos%z  .and. 
     *         Volat( Det%cmp(ncx)%vol + pipeh) .ge. pos%z) then
          temp = pos%x**2 + pos%y**2
          if(temp .le. Volat( Det%cmp(ncx)%vol + pipeor)**2 .and.
     *       temp .ge. Volat( Det%cmp(ncx)%vol + pipeir)**2)  then
             icon = 0
          else
             icon = 1
          endif
       else
          icon = 1
       endif
       end
       subroutine epssbox(ncx, pos, icon)
!          search specified box
!          Is pos in ncx-th comp.  ? if yes, icon=0 else
!          icon = 1
!          pos is assumed to be in Local coord.
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos
           integer ncx,  icon
        

           if(0.d0 .le. pos%z  .and. 
     *        pos%z .le. Volat( Det%cmp(ncx)%vol + boxc)) then
              if(0.d0 .le. pos%x .and.
     *             pos%x .le. Volat( Det%cmp(ncx)%vol + boxa)) then
                 if(0.d0 .le. pos%y .and.
     *                pos%y .le. Volat( Det%cmp(ncx)%vol + boxb)) then
                    icon  =0
                    goto 101
                 endif
              endif
           endif
           icon = 1
 101       continue
       end
       subroutine epsssph(ncx, pos, icon)
!          search specified sphere
!         Is  pos in ncx-th comp.  ? if yes, icon=0 else
!               icon = 1
!         pos is assumed to be Local coord.
       implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos
           integer ncx,  icon

           if( pos%x**2 + pos%y**2 + pos%z**2  .le.
     *         Volat( Det%cmp(ncx)%vol + sphr)**2) then
              icon  = 0
           else
              icon = 1
           endif
       end
       subroutine epw2l(ncx, pos, poso)
!        world to local conversion
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"

       type(epPos)::  pos    ! input
       type(epPos)::  poso    ! output. poso can be pos
       type(epPos)::  temp
      integer ncx
      if(ncx > Det%nct) then
         poso = pos
         return !**********
      endif

      poso%x = pos%x - Det%cmp(ncx)%orgx
      poso%y = pos%y - Det%cmp(ncx)%orgy
      poso%z = pos%z - Det%cmp(ncx)%orgz
 
      if(Det%cmp(ncx)%rotation) then
         temp%x = Det%cmp(ncx)%direc(1) * poso%x +
     *      Det%cmp(ncx)%direc(2) * poso%y +
     *      Det%cmp(ncx)%direc(3) * poso%z


         temp%y = Det%cmp(ncx)%direc(4) * poso%x +
     *      Det%cmp(ncx)%direc(5) * poso%y +
     *      Det%cmp(ncx)%direc(6) * poso%z

         temp%z = Det%cmp(ncx)%direc(7) * poso%x +
     *      Det%cmp(ncx)%direc(8) * poso%y +
     *      Det%cmp(ncx)%direc(9) * poso%z
         poso = temp
      endif

      end
      subroutine epl2w(ncx, pos, poso)
!         local to world conversion
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epPos)::  pos ! input
       type(epPos)::  poso ! output, poso can be pos

      integer ncx
       type(epPos)::  temp

      if(Det%cmp(ncx)%rotation) then
         temp%x = Det%cmp(ncx)%direc(1) * pos%x +
     *      Det%cmp(ncx)%direc(4) * pos%y +
     *      Det%cmp(ncx)%direc(7) * pos%z 


         temp%y = Det%cmp(ncx)%direc(2) * pos%x +
     *      Det%cmp(ncx)%direc(5) * pos%y +
     *      Det%cmp(ncx)%direc(8) * pos%z

         temp%z = Det%cmp(ncx)%direc(3) * pos%x +
     *      Det%cmp(ncx)%direc(6) * pos%y +
     *      Det%cmp(ncx)%direc(9) * pos%z

         poso%x = temp%x + Det%cmp(ncx)%orgx
         poso%y = temp%y + Det%cmp(ncx)%orgy
         poso%z = temp%z + Det%cmp(ncx)%orgz
      else
         poso%x = pos%x + Det%cmp(ncx)%orgx
         poso%y = pos%y + Det%cmp(ncx)%orgy
         poso%z = pos%z + Det%cmp(ncx)%orgz
      endif


      end
      subroutine epw2ld(ncx, pos, poso)
!        world to local conversion for direction cos.
!        also you can use this vector direction such as B, E
!        field.
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epDirec)::  pos ! input
       type(epDirec)::  poso ! output. poso can be pos

      integer ncx
 
       type(epDirec)::  temp
!/////////
!      write(0,*) ' in epw2ld . ncx=',ncx, ' Det.nct=',Det.nct
!      write(0,*) ' pos=',pos.x,pos.y, pos.z
!/////////////
      if(ncx > Det%nct) then
         poso = pos
         return   !*********
      endif
      if(Det%cmp(ncx)%rotation) then
         temp%x = Det%cmp(ncx)%direc(1) * pos%x +
     *      Det%cmp(ncx)%direc(2) * pos%y +
     *      Det%cmp(ncx)%direc(3) * pos%z


         temp%y = Det%cmp(ncx)%direc(4) * pos%x +
     *      Det%cmp(ncx)%direc(5) * pos%y +
     *      Det%cmp(ncx)%direc(6) * pos%z

         temp%z = Det%cmp(ncx)%direc(7) * pos%x +
     *      Det%cmp(ncx)%direc(8) * pos%y +
     *      Det%cmp(ncx)%direc(9) * pos%z

         poso = temp
      else
         poso = pos
      endif


      end
      subroutine epw2ldm(ncx, pos, poso, pinout)
!        world to local conversion for direction cos and
!        momentum
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epDirec)::  pos ! input
       type(epDirec)::  poso ! output. poso can be pos
       type(ptcl)::  pinout ! in/out ptcl info.
      integer ncx
      real*8 pabs
       type(epDirec)::  temp

      if( ncx > Det%nct ) then
         poso = pos
         return  !*********
      endif
      if(Det%cmp(ncx)%rotation) then
         temp%x = Det%cmp(ncx)%direc(1) * pos%x +
     *      Det%cmp(ncx)%direc(2) * pos%y +
     *      Det%cmp(ncx)%direc(3) * pos%z


         temp%y = Det%cmp(ncx)%direc(4) * pos%x +
     *      Det%cmp(ncx)%direc(5) * pos%y +
     *      Det%cmp(ncx)%direc(6) * pos%z

         temp%z = Det%cmp(ncx)%direc(7) * pos%x +
     *      Det%cmp(ncx)%direc(8) * pos%y +
     *      Det%cmp(ncx)%direc(9) * pos%z
         poso = temp
         pabs = sqrt(pinout%fm%p(1)**2 +pinout%fm%p(2)**2
     *          + pinout%fm%p(3)**2)
         pinout%fm%p(1) = pabs*temp%x
         pinout%fm%p(2) = pabs*temp%y
         pinout%fm%p(3) = pabs*temp%z
      else
         poso = pos
      endif


      end
      subroutine epl2wd(ncx, pos, poso)
!         local to world conversion for direction cos
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epDirec)::  pos ! input
       type(epDirec)::  poso ! output. poso can be pos

      integer ncx
       type(epDirec)::  temp

      if(Det%cmp(ncx)%rotation) then
         temp%x = Det%cmp(ncx)%direc(1) * pos%x +
     *      Det%cmp(ncx)%direc(4) * pos%y +
     *      Det%cmp(ncx)%direc(7) * pos%z 


         temp%y = Det%cmp(ncx)%direc(2) * pos%x +
     *      Det%cmp(ncx)%direc(5) * pos%y +
     *      Det%cmp(ncx)%direc(8) * pos%z

         temp%z = Det%cmp(ncx)%direc(3) * pos%x +
     *      Det%cmp(ncx)%direc(6) * pos%y +
     *      Det%cmp(ncx)%direc(9) * pos%z
         poso = temp
      else
         poso = pos
      endif


      end
      subroutine epl2wdm(ncx, pos, poso, pinout)
!         local to world conversion for direction cos
!          and momentum
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
       type(epDirec)::  pos ! input
       type(epDirec)::  poso ! output. poso can be pos
       type(ptcl)::   pinout !  for momentum part are affected

      integer ncx
       type(epDirec)::  temp
      real*8 pabs

      if(Det%cmp(ncx)%rotation) then
         temp%x = Det%cmp(ncx)%direc(1) * pos%x +
     *      Det%cmp(ncx)%direc(4) * pos%y +
     *      Det%cmp(ncx)%direc(7) * pos%z 

         temp%y = Det%cmp(ncx)%direc(2) * pos%x +
     *      Det%cmp(ncx)%direc(5) * pos%y +
     *      Det%cmp(ncx)%direc(8) * pos%z

         temp%z = Det%cmp(ncx)%direc(3) * pos%x +
     *      Det%cmp(ncx)%direc(6) * pos%y +
     *      Det%cmp(ncx)%direc(9) * pos%z
         poso = temp
         pabs = sqrt( pinout%fm%p(1)**2 + pinout%fm%p(2)**2
     *               + pinout%fm%p(3)**2 )
         pinout%fm%p(1) = pabs*temp%x
         pinout%fm%p(2) = pabs*temp%y
         pinout%fm%p(3) = pabs*temp%z
      else
         poso = pos
      endif
      end

      subroutine epv2c_cyl(comp, posv, posc)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  posv ! input
       type(epPos)::  posc ! output. 
       type(Component)::  comp
      if( comp%struc(1:5) == "cyl_x") then
!         posc = epPos( posv.z, posv.y, posv.x)
         posc = epPos( posv%y, posv%z, posv%x)
      elseif(comp%struc(1:5) == "cyl_y") then
!         posc = epPos( posv.x, posv.z, posv.y)
         posc = epPos( posv%z, posv%x, posv%y)
      else
         posc = posv
      endif
      end
      subroutine epc2v_cyl(comp, posc, posv)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  posc ! input
       type(epPos)::  posv ! output. 
       type(Component)::  comp
      if( comp%struc(1:5) == "cyl_x") then
!         posv = epPos( posc.z, posc.y, posc.x)
         posv = epPos( posc%z, posc%x, posc%y)
      elseif(comp%struc(1:5) == "cyl_y") then
!         posv = epPos( posc.x, posc.z, posc.y)
         posv = epPos( posc%y, posc%z, posc%x)
      else
         posc = posv
      endif
      end
      subroutine epv2cd_cyl(comp, dirv, dirc)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  dirv ! input
       type(epPos)::  dirc ! output. 
       type(Component)::  comp
      if( comp%struc(1:5) == "cyl_x") then
!         dirc = epPos( dirv.z, dirv.y, dirv.x)
         dirc = epPos( dirv%y, dirv%z, dirv%x)
      elseif(comp%struc(1:5) == "cyl_y") then
!         dirc = epPos( dirv.x, dirv.z, dirv.y)
         dirc = epPos( dirv%z, dirv%x, dirv%y)
      else
         dirc = dirv
      endif
      end
      subroutine epc2vd_cyl(comp, dirc, dirv)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  dirc ! input
       type(epPos)::  dirv ! output. 
       type(Component)::  comp
      if( comp%struc(1:5) == "cyl_x") then
!         dirv = epPos( dirc.z, dirc.y, dirc.x)
         dirv = epPos( dirc%z, dirc%x, dirc%y)
      elseif(comp%struc(1:5) == "cyl_y") then
!         dirv = epPos( dirc.x, dirc.z, dirc.y)
         dirv = epPos( dirc%y, dirc%z, dirc%x)
      else
         dirv = dirc
      endif
      end

      subroutine epv2c_pipe(comp, posv, posc)
!           variant local to canonical local; used when
!           manipulating things in canonical local
!          (for boundary search etc).  system normally recognizes
!           variant local as usual local 
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  posv ! input
       type(epPos)::  posc ! output. 
       type(Component)::  comp
      if( comp%struc == "pipe_x") then
!         posc = epPos( posv.z, posv.y, posv.x)
         posc = epPos( posv%y, posv%z, posv%x)
      elseif(comp%struc == "pipe_y") then
!         posc = epPos( posv.x, posv.z, posv.y)
         posc = epPos( posv%z, posv%x, posv%y)
      else
         posc = posv
      endif
      end
      subroutine epc2v_pipe(comp, posc, posv)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  posc ! input
       type(epPos)::  posv ! output. 
       type(Component)::  comp
      if( comp%struc == "pipe_x") then
!         posv = epPos( posc.z, posc.y, posc.x)
         posv = epPos( posc%z, posc%x, posc%y)
      elseif(comp%struc == "pipe_y") then
!         posv = epPos( posc.x, posc.z, posc.y)
         posv = epPos( posc%y, posc%z, posc%x)
      else
         posc = posv
      endif
      end
      subroutine epv2cd_pipe(comp, dirv, dirc)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  dirv ! input
       type(epPos)::  dirc ! output. 
       type(Component)::  comp
      if( comp%struc == "pipe_x") then
!         dirc = epPos( dirv.z, dirv.y, dirv.x)
         dirc = epPos( dirv%y, dirv%z, dirv%x)
      elseif(comp%struc == "pipe_y") then
!         dirc = epPos( dirv.x, dirv.z, dirv.y)
         dirc = epPos( dirv%z, dirv%x, dirv%y)
      else
         dirc = dirv
      endif
      end
      subroutine epc2vd_pipe(comp, dirc, dirv)
      implicit none
#include "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"
       type(epPos)::  dirc ! input
       type(epPos)::  dirv ! output. 
       type(Component)::  comp
      if( comp%struc == "pipe_x") then
!         dirv = epPos( dirc.z, dirc.y, dirc.x)
         dirv = epPos( dirc%z, dirc%x, dirc%y)
      elseif(comp%struc == "pipe_y") then
!         dirv = epPos( dirc.x, dirc.z, dirc.y)
         dirv = epPos( dirc%y, dirc%z, dirc%x)
      else
         dirv = dirc
      endif
      end
      subroutine epGetBaseStrucName(struc, mname)
      implicit none
      character(*),intent(in)::struc  ! box, box_w cyl cyl_y_w etc
      character(*),intent(out):: mname  ! box cyl etc without _
!            normally   their length = MAX_STRUCCHR but could be less or more
!        
      integer:: i
      i = index(struc, '_')
      if( i > 0 ) then
         mname=struc(1:i-1)
      else
         mname= struc
      endif
      end   subroutine epGetBaseStrucName
