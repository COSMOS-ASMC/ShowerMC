!      subroutine epbndry(pos, icon)
      subroutine epbndry(pos, el, icon)
      use moddebug
      implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
!
!       current track pos is assumed to be inside the current
!       component. 
!    At first, ptcl is inside the world. Then, this routine tries
!   to find x-point with the wrold itself. It should exit. 
!   However, there may be  contained components. So x-points
!   with all ingredients are searched for by using epbndry2
!   placed not at the first part of epbndry. The nearest 
!   x-point is selected and return is made with icon = 0.
!   The ptcl point is moved a little bit (EpsLeng) from the
!   X-ing point, so it must be inside of the new component.
!
!       get the crossing point of cTrack extending from the current
!       location with a given direction
!
       type(epPos)::  pos ! output. obtained crossing point 
                         ! in the current local coord.
      integer icon   !  output 0 --> normal. pos, el found
                     !         1 --> something wrong. adjuested safely.
                     !               pos. el obtained  
                     !         2---> ptcl is near void and judged
                     !               inside at some place and
                     !               in void at other place. so 
                     !               discard it.
      character*120 msg
      real*8  tmpel, el2, el
      integer i, cnxx, j, cnx, jcon
       type(epPos):: postemp
       type(epDirec)::  dirtemp
       type(epTrack)::  saveTrack
      integer,save:: PrevNo
      save saveTrack
      character(len=MAX_STRUCCHR):: basename
!          present comp, boundary
      call epbndry2(Cn, el, jcon)
      
      if(jcon .ne. 0) then
!               No boundary found in the  present comp. strange.
!         This happens if a partcle is on a boundary or 
!         very near to the boudnary.  We can adust the  position
!         etc and continue the execution.
!         However, this also happens if a wrong incident direcion is given,
!         too or by some unkown reason. So we count the occurence / event
!         issue messages.
         if( Bndryerr == 0 ) then
            ! save pos and dir
            saveTrack = cTrack
            PrevNo = 0
         endif
         Bndryerr = Bndryerr + 1    ! count it
         PrevNo = PrevNo + 1
         call eppos2cn(Cn, cTrack, cnxx)
         if(cnxx .ge. Det%nct) then
            icon = 2
            return !  *******
         endif
         if( Bndryerr < 10 ) then
            call epGetBaseStrucName(Det%cmp(Cn)%struc, basename)
!              sqtccl has some problem and when it is treated
!              in some case error happens; it is difficult to
!              fix completely due to numerical error. It is
!              proved to be non harmul (atuomatically remedied)
!              so avoid message  until complete fix comes.
!              But see later, if it is really problem, 
!              next call will be done for sqtccl too.
            if( basename /= "sqtccl") then
               call epIssueBndryErrMsg
!!!!!!!!!!!!!!!!!
!　　　　　　   call cerrorMsg('hcs+boundary', 0)
!               stop                 
!!!!!!!!!!!!               
            endif
         endif

         if( PrevNo  > 1 ) then
            if( cTrack%p%code /= saveTrack%p%code .or.
     *          cTrack%p%charge /= saveTrack%p%charge .or.
     *          cTrack%p%subcode /= saveTrack%p%subcode ) then
                prevNo = 0
                saveTrack = cTrack
             elseif( 
     *          abs(cTrack%pos%x-saveTrack%pos%x) > EpsLeng .and.
     *          abs(cTrack%pos%y-saveTrack%pos%y) > EpsLeng .and.
     *          abs(cTrack%pos%z-saveTrack%pos%z) > EpsLeng  ) then
                PrevNo = 0
                saveTrack = cTrack
             endif
          endif
          if( PrevNo > 10 ) then
              ! hopeless to escape err; discard this
             if( cTrack%p%fm%p(4) < Total1ryE* 1.d-4) then
                icon = 2
                return          !  *******
             endif
          endif
             
         if( Bndryerr  > 
     *     min(500,  max(int(10*Total1ryE/100), 10)) ) then
           if(basename == "sqtccl")  call  epIssueBndryErrMsg
            call cerrorMsg('too many boundary error', 0)
         endif
         if(Cn .ne. cnxx) then
            cTrack%cn = cnxx
            call epl2w(Cn, cTrack%pos, postemp)
            call epw2l(cnxx, postemp,  cTrack%pos)

!!          call epl2wd(Cn, cTrack%w,  dirtemp)   !   2016/Sep 
!1                      replaced by next
            call epl2wdm(Cn, cTrack%w,  dirtemp, cTrack%p)
            
            call epw2ldm(cnxx, dirtemp, cTrack%w, cTrack%p)
            pos = cTrack%pos
            call epnewComp(cTrack)
         else
            pos = cTrack%pos
              ! go back little bit
!            cTrack.pos.x = cTrack.pos.x - EpsLeng*cTrack.w.x
!            cTrack.pos.y = cTrack.pos.y - EpsLeng*cTrack.w.y
!            cTrack.pos.z = cTrack.pos.z - EpsLeng*cTrack.w.z
            dirtemp%x=1.
            dirtemp%y=0.
            dirtemp%z = 0.
            call ctransVectZ(cTrack%w, dirtemp, dirtemp)
            cTrack%pos%x = cTrack%pos%x - EpsLeng*dirtemp%x
            cTrack%pos%y = cTrack%pos%y - EpsLeng*dirtemp%y
            cTrack%pos%z = cTrack%pos%z - EpsLeng*dirtemp%z
         endif
         el = EpsLeng
         icon = 1
         jcon = 0
         CrossMode = 0
      else
         CrossMode = 0
!         **** partial**** container's boundary,
!               if any, may be closer
         do j = 1, Det%cmp(Cn)%NPContainer
!((((((((((
!           cnx = Det.cmp(Cn).PContained(j)
            cnx = CnArea(Det%cmp(Cn)%PContained+j)
!))))))))))))
!!!!!!!///////
!            write(0,*) ' calling bndry2 B'
!////////////
            call  epbndry2(cnx,  el2, jcon)
            if(jcon  .eq. 0) then
               if(el2 .lt.  el) then
                  el = el2
                  CrossMode = 2
               endif
            endif
         enddo
!         *** Matreska *** may be in between
         do i = 1, Det%cmp(Cn)%NMatreska
!((((((((((((((((
!!!!!!!///////////
            call epbndry2(CnArea( Det%cmp(Cn)%Contains+i), tmpel, jcon)
!))))))))))
            if(jcon .eq. 1 .and. tmpel .lt. el) then
               el = tmpel
               CrossMode = 1
            endif
         enddo
         pos%x = cTrack%pos%x + el*cTrack%w%x
         pos%y = cTrack%pos%y + el*cTrack%w%y
         pos%z = cTrack%pos%z + el*cTrack%w%z
         icon = 0
      endif
      end

      subroutine epIssueBndryErrMsg
      implicit none
      call cerrorMsg(' ',  1)
      call cerrorMsg(
     *  'If following message comes out at the very beginning of', 1)
      call cerrorMsg(
     *  ' an event, better to suspect that your incident particle',1)
      call cerrorMsg(
     *  'direction might be wrong; check it.!!!!!', 1)
      call cerrorMsg(
     *  'If not,  a particle on or very near to a ', 1)
      call cerrorMsg(
     *  'component boundary might have buffalloed Epics.' , 1) 
      call cerrorMsg(
     *  'Epics is expected to recover the matter safely.', 1)
      call cerrorMsg(
     *  'However, some unknown cause may lead to this, ', 1)
      call cerrorMsg(
     *  'so check the coordinate below, if they are very ',1)
      call cerrorMsg(
     *  'near to the boundary, it will  be ok', 1)
      call epfordebug('epbndry')
      call cerrorMsg(
     *   'NOTE: local above means local coordinate', 1)
      call cerrorMsg(
     *   'so if component name has "_", it is different', 1)
      call cerrorMsg( '  from  canonical coordinate', 1)
      call cerrorMsg('-------------------------------',1)
      end   subroutine epIssueBndryErrMsg

