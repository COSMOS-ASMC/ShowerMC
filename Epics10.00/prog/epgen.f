!     those moved to Package/EM/FromEpi are still exist here so remove them later
!      
#include "ZsubstRec.h"
#include "Zepcondc.h"
      module moddebug
      logical,save:: debugmodeON=.false.
      logical,save:: firstdebug=.true.
      end module moddebug
!        epifCross is changed radically.
!        epbndry; icon is changed to get el in real*8
!        *********************************************************
!        *
!        *  epgen:generate  showers until stack area becomes empty
!        *
!        *********************************************************
!
      module moddedx
!         work common block for dedx (GeV/g/cm2) (epEloss / epchcke)
      real(8):: dedx  !  restricted dE/dx 
      real(8):: dedxf  ! full dE/dx
      end module moddedx

      subroutine epgen
      use modV1ry
      use modMCScontrol
      implicit none
#include  "Zmedia.h"
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcnfig.h"
#include  "ZsepManager.h"
! >>>>>>>>>>>>>>>light
      integer,save::onlyonce
!      common /dddebug/ onlyonce
!<<<<<<<<<<<<<
      integer icon
!         init. for 1 event has been finished, next is
!       to put other final init. for 1 event

      call epr1ev
!>>>>>>>>>>>>>>>>>>>light
      onlyonce = 0
!<<<<<<<<<<<<<<<<<<
!         *** until loop*** 
      icon = 0
      do while (icon .eq. 0)
!                get 1 particle from stack area 
         call eppop(cTrack, icon)
!                icon=0: 1 ptcl gotten in cTrack
!                       -1: no more ptcl  in the stack
         if(icon .eq. 0) then
            if( V1ry /= 0 ) then
               if( cTrack%inciflag  == 1 ) then
                  V1ry = abs(V1ry)
               else
                  V1ry = -abs(V1ry)
               endif
            endif
            call epfl1
         else
           !  all of normal tracking finished.>>>>>>>>>>>>>>>>light
            if( Light == 12 .and. onlyonce == 0 ) then
              ! energy deposit has been stored in cells.  we must generate
              ! scinti light and do ray tracing.  We first push Edepo as
              ! psudo ptcls
               call epLightPushCells(0, 0)
               icon = 0  ! now many cells stacked; reset icon
               onlyonce = 1
            endif   
           !<<<<<<<<<<<<<<<<<<<<<<<<<<
         endif
      enddo
       !>>>>>>>>>>>>>>>>>>>>>>>light
      if( Light == 21 ) then
         ! cell stored deposit should be put in a disk  as primary
         ! for another job (file name with + or -)
         if( FirstC ) then
            ! no interatcion of incident so far; we must write
            ! event header
            call epLightIOwrite1stCol
         endif
         call epLightIOwriteCell
         call epLightIOwritedE(OutPrimEff, Det%nct)
             ! enery deposit of all comp. with d /=0
      endif
       !<<<<<<<<<<<<<<<<<<<<<<
      end
      subroutine epicosmos(param)
      implicit none
!          read parameter file for Cosmos.
#include  "Zmanagerp.h"
      integer  icon
      character*(*) param       ! input. Cosmos param file path.
      call copenNLf(TempDev, param, icon)
      if(icon .ne. 0) then
         call cerrorMsg('epicosmos cannot open parameter file',0)
      endif
      call creadParam(TempDev)
      close(TempDev)
      end
      
      subroutine epiaev(dsn1, dsn2)
      implicit none
!             init for all event; read epics data and cofig data
#include  "ZepManager.h"
#include  "ZsepManager.h"
#include  "Zmedia.h"      
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcnfig.h"


      character*(*) dsn1        ! input epics data file path
      character*(*) dsn2        ! input config data file path

      integer klena
      character*16 uid
      integer i, icon
!
!                 open basic data residence file for epics
!              component # is undefined yet
      Cn=-1
      call epempty
!      Stack_pos = 0
!      StackDisk = 0
      Nevrun=1   !  >= 9.06 (<=9.05; Nevrun=0)
      DoUpdateNo =.true.  ! >=9.06
!                read basic parameter

      call epprmr(dsn1)
!                compute something using basic parameters
      call epcmp1

!                read configuration

!     
       call eprcnf(dsn2)
       call epiniXsMedia  ! make media info needed for H.I, MCS
                         ! under Cosmos
       call cHowMCS  ! MCSmodel is analysed ...
       call ciniMCS  ! init for MCS and set doNewMCS (t/f)


!      *********************
!               init for nuclin, hadrin
!c       call haddenC
!c       call chanwnC
         call cintModels('epics')    ! this is probalby not needed now
!       ################
!                user dependent init
       call uiaev
!>>>>>>>>>>>>>>>>light
!         if Light >0, manipulate CountDE of each component
       if(Light >  0) then
          call epLightCountDE ! 0 means from system

          call epLightAlloc  ! alloc arrays
          if(Light == 21) then
             call epLightIOwriteIni
          endif
       endif
!<<<<<<<<<<<<<<
       !  allocate Eloss array for all comp. (always)
       call epAlloc(Det%nct)

       if(MagField .eq. 1)then
          Bfield%x = Bxu
          Bfield%y = Byu
          Bfield%z = Bzu
       endif
       if(ElecField .eq. 1) then
          Efield%x = Exu
          Efield%y = Eyu
          Efield%z = Ezu
       endif
!
!          fix tracedir ; This must be before uiaev in future
       if(Trace) then
          if(TraceDir .eq. ' ') then
             call cgetLoginN(uid) ! cosmos function
             TraceDir ='/tmp/' // trim(uid)
          endif
       endif
       end
       subroutine  epi1ev(icon)
!             init for 1 event
       use modIntInfo
       implicit  none
#include "ZepManager.h"
#include "ZepTrackp.h"
#include  "Zmedia.h"       
#include "ZepTrackv.h"

       integer icon ! output. if  0,  user has set ptcl in stack
                   !             1,  no ptcl has set in stack
      character*200 msg
      character*128 filen
      character*16 uid
      integer jcon, leng
      integer  klena, i

       Bndryerr = 0   ! counter for bundary search failures/ event
       Making1ry = 1  ! Now preparing 1 event so 1ry is being made
                    ! will be made to be 0 at epr1ev.
                    ! If this is 1, partcles being stacked are 1ry
                    ! we counts total energy of 1ry.
       Total1ryE = 0. ! clear total 1ry K.E counter of the current event
       Cn=-1
       Firsti=.true.
       FirstC=.true.
       Proc1 = '   '   ! first collision process.
       FirstInt%x = -100000.d0
       FirstInt%y = -100000.d0
       FirstInt%z = -100000.d0
       FirstCn = 0 
       FirstMedia%colElem = 0

#if defined (INTINFO)
!       By  next, the user interface epUI is called 
       codeAforInt(:)= 0
#endif


!       Nevrun = Nevrun + 1     ! <= 9.05
       pathInB=0.
       SumDe = 0.
       Move%Abort = 0     !  if this become non zero for an event, 
                     !  the particle is discarded or event is discarded
       call epcEloss  ! clear Eloss counters
!>>>>>>>>>>>>>>>Light
       if(Light > 0) then
          call epLighti1ev
          if(Light == 21 ) then
             call epLightIOwriteiev
          endif
       endif
!<<<<<<<<<<<<<<<<<
!                 user dependent init
       call ui1ev
!               next call sets to icon=1,    if usetip is not
!               user made.
       call usetip(icon)
!                 icon=0-->user routine already set prtcl
!                 icon=1-->no user routine for incident
       if(Trace) then
          write(filen, *) trim(TraceDir)//'/trace',
     *         Nevrun 
          call kseblk(filen, ' ', leng)
          call copenfw(abs(IoTrace), filen,  jcon)
          if(jcon .ne. 0) then
             write(0,*)
             call cerrorMsg(
     *       'You gave Trace=t in sepicsfile, but the file ', 1)
             write(0,*) trim(filen)
             write(msg, '(a,a,a)')
     *       ' cannot be opened: Probably you have to make', 
     *         trim(uid), ' directory. Or '
             call cerrorMsg(msg, 1)
             write(0,*) 'If you gave TraceDir in both "epicsfile" and' 
             write(0,*) '"sepicsfile", the former overrides the latter'
             write(0,*) 'so the problem may be solved by removing ' 
             write(0,*) 'TraceDir in epicsfile'
             stop 
          endif
       endif

       end
!      *******************
       subroutine  epe1ev
!      ******************* end of 1 event
       implicit none
#include "ZepManager.h"
#include  "ZepTrackp.h"
#include  "Zmedia.h"       
#include "ZepTrackv.h"

!>>>>>>>>>>>>>>>>>Light
       if(Light >  0) then
          call epLighte1ev
       endif
!<<<<<<<<<<<<<<<<<<<
!               user dependent end process

       if(Move%Abort .le. 1) call ue1ev
       if( DoUpdateNo ) then 
          Nevrun = Nevrun +1
       else
          DoUpdateNo = .true.
       endif
       if(Trace) then
          close(abs(IoTrace))
       endif
       end
!      *************  now 1 event is ready to start
       subroutine epr1ev
       implicit none
#include  "Zmedia.h"       
#include "ZepTrackv.h"
#include "ZepTrackp.h"

       integer icon

!          the procedures to be performed after the completion of
!          init. for 1 event must be placed here
!          save first stacked track as  incident
       if( Light /= 22 ) then
               ! in the case of 22, primary info has been
               ! read alredy, since it is placed at the top
               ! of each event
          call epgetTrack(1, Incident, icon)
       endif
!        if(icon .ne. 0) then
!           call cerrorMsg('no incident found in epr1ev', 0)
!        endif
!
        call uafi1ev            ! after init of 1 event
        Making1ry = 0        
        end
!      ****************************
!                user dependent all event end
       subroutine epeaev
       implicit none

       call ueaev
       end

!c      ***********************
      subroutine epSkipUpdateNo
      implicit none
#include  "ZepManager.h"
!            disable the update of Nevrun
      DoUpdateNo = .false.

      end

      subroutine epfl1
!!!!!
      use moddebug
!!!!!!!!!      
       implicit none
!           follow 1 particle
!              until current partile dies or all are put into stack.
#include  "ZepTrackp.h"
#include  "Zmedia.h"       
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"

      integer icon, info

      do while (.true.)
!                current track is in  cTrack
!          fix process and free path
!     and compute tentative new pos.--> Move
         call epnewp(icon)
         if(icon .ne. 0) goto 100 ! abort is specified by the user
!          energy loss count in user hook
!         if(Det.cmp(Cn).CountDE .gt. 0 .and. >>>>>>>>>light
          if(  cTrack%p%charge .ne. 0) then  ! <<<<<<<<< 
!             if( Move.Track.p.fm.p(4)-Move.Track.p.mass
!     *                                     .le. KEmin ) then
!                info = 1
!             else
                info = 0
!             endif
!               info may be used only to set Move.Abort    
             call epLightPreUserde(info, cTrack)
             if(Move%Abort /= 0 ) then
                if(Move%Abort /= 3 ) then
                   call  epempty ! empty the sack. discard ev.
                   call  epSkipUpdateNo
                else
                   Move%Abort = 0
                endif
                icon = 1
                goto 100
             endif

          endif
!          add  time  
         if(TimeStruc) then
            call epaddTime
         endif

!           adjust momentum; because of  energy change
         if(cTrack%p%charge .ne. 0) then
            call epe2p(Move%Track)  ! not Move.Track.p; bug <=v9.13
         endif

!            if(.not. Move.Cross .and. cTrack.p.charge .ne. 0) then
!               we don't use above judgement; later at epCross, we use new
!               angle due to scattering at Move.boundary. (don't use
!               new  position due to scattering in epCross)
!           ++++++++++
         if(cTrack%p%charge .ne. 0) then

!                 multiple scattering, magnetic deflection
!                 electric deflection. 
!     if MCSmode, MCS will not be done inside
            call epdeflection(icon)  
         endif
!              take trace info.
         call epTraceFE  ! v9.17
!         if(Trace) then
!            if(IoTrace .lt. 0 .or. (IoTrace .gt. 0 
!     *           .and. cTrack.p.charge .ne. 0)) then
!               call epTrace
!            endif
!         endif
!              really move
#ifdef  SUBSTREC
         cTrack = Move%Track
#else
         call epsubstTRK( cTrack, Move%Track)
#endif
         if(.not. Move%Trunc) then
            call epint(icon)    ! interaction routine
!                icon =  1 always.
         elseif(Move%Cross) then
            call epCross(icon)
         endif
         if         (icon .ne. 0)
     *                       goto 100
      enddo
 100  continue
      end

!     *****************
      subroutine epnewp(icon)
!     *****************
!   1)    copy current track to Move.track
!   2)    compute lenghtToB= length (cm) to the boundary
!          or crossing point with contained component
!          if somethig very strange, icon =2 and return
!   3)    chckE0: simple energy check. 
!      If E< Emin, 
!         3-1) compute the ragne of the ptcl.
!             if the range < lengthToB, absorbe energy
!                   and return with icon = 1
!             elseif( Ek < AEmin ) , abosrbe energy 
!                   and return with icon =1  ; to avoid
!                   delicate problem when lengthToB ~0 and E is
!                   very small
!             endif
!             do as if E> Emin (go next) 
!      endif
!
!                
       use modepRange, only:       attenF
       use modMCScontrol
       use modV1ry
!!!!!!!!
       use moddebug
!!!!!!!!       
      implicit none
#include  "ZepManager.h"
#include  "ZmediaLoft.h"      
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"
      
      integer icon
      logical Absorb
      real(8)::AbsoEmin=2.d-6
      real(8):: Rgm, Rcm
      real(8):: u, fpos

!    *********************** sample interaction length ****************

!        copy current track  to Move
#ifdef  SUBSTREC
      Move%Track = cTrack
#else
      call epsubstTRK( Move%Track, cTrack)
#endif
!        get the crossing point of the current track
!        with the current component or contained component
!        el is the length to the crossing point  
      call epbndry(Move%boundary, lengthToB, icon)
      if(icon == 2 )  return  ! *****
!       icon =1 indicates something wrong. but 
!       it should have been corrected safely so we don't
!       care.

      call epchckE0(cTrack, icon)  ! simple energy check
      if( icon /= 0 ) then
         if( ImperativeEmin .or. icon == -1 ) then
            Absorb = .true.
            Move%dl = 0.
            Move%dx = 0.
            Move%dt = 0.
         elseif(( AutoEmin == 2 .or. AutoEmin == 4 ) .or.
     *     (cTrack%p%code == kphoton .and. !  for photon v9.154
     *       AutoEmin == 3 ) ) then
!             energy is low, see if the range is < el
            call epGetRange(MediaNo, Media(MediaNo), cTrack%p,
     *                      Rgm, Rcm)
            Absorb =  Rcm < lengthToB
            if( Absorb ) then
            !  current ptcl KE may be absorbed during path = Rcm
            !  suppose particle can run Rcm (for charge).
            !  For photons, sample the absroption point by
            !  using attentuation length (=Rcm/attenF)
               if( cTrack%p%code == kphoton ) then
                  call rndc(u)
                  fpos=min( -log(u)/attenF, 1.d0)  ! l
                  Rcm = Rcm * fpos
                  Rgm = Rgm * fpos 
               endif
               Move%dl = Rcm
               Move%dx = Rgm
               Move%dt = Move%dx/Media(MediaNo)%X0g
            elseif ( cTrack%p%fm%p(4) - cTrack%p%mass < AbsoEmin) then
              ! for e, very rare to come.  neutrons are more 
              ! others are very very rare
               Absorb = .true.
               Move%dl = lengthToB/3.
               Move%dx = Move%dl/Media(MediaNo)%gtocm
               Move%dt = Move%dx/Media(MediaNo)%X0g
            endif
         else
            Absorb = .true.
            Move%dl = 0.
            Move%dx = 0.
            Move%dt = 0.
         endif
         if( Absorb ) then
            Move%Track%p%fm%p(4) = Move%Track%p%mass
            Move%Track%pos%x = cTrack%pos%x + cTrack%w%x*Move%dl
            Move%Track%pos%y = cTrack%pos%y + cTrack%w%y*Move%dl
            Move%Track%pos%z = cTrack%pos%z + cTrack%w%z*Move%dl
            call epAbsorb(cTrack, icon)  ! icon = 1
            call epTraceFE    ! v9.17 
!            if(Trace) then
!               if(IoTrace .lt. 0 .or. (IoTrace .gt. 0
!     *              .and. cTrack.p.charge .ne. 0)) then
!!                      next moved above
!!                  Move.Track.p.fm.p(4) = Move.Track.p.mass
!!                  Move.Track.pos.x = cTrack.pos.x + cTrack.w.x*Move.dl
!!                  Move.Track.pos.y = cTrack.pos.y + cTrack.w.y*Move.dl
!!                  Move.Track.pos.z = cTrack.pos.z + cTrack.w.y*Move.dl
!                  call epTrace
!               endif
!            endif
            icon = 1
            return
         endif
      endif
!           ptcl still alive;  interaction path     
!       do for each ptcl

!      ActiveMCS = 'Mol' ==>
      call csetActiveMCS
      
      call ciniIntInf           !    clear process recorder
      
      if(cTrack%p%code .eq. kphoton) then
         call epprog(Media(MediaNo))  ! Move.proc, Move.dt fixed         
      elseif(cTrack%p%code .eq. kelec) then
         call epproe(Media(MediaNo))            ! Move.proc, Move.dt fixed
      elseif( cTrack%p%code > 0 ) then
         call epNonEleMag(Media(MediaNo))
!>>>>>>>>>>>>>light
      elseif( Light > 0 ) then
         if( cTrack%p%code == klight ) then
            call epproLight
         else
            write(0,*) " invalid  code=",cTrack%p%code
            write(0,*) " for light related ptcl's path sampling"
            stop
         endif
      endif
!     *************** truncate if too long  **********************
!              trancate if dt is too long (set Trunc=t/f)
      call eptrunc
      if( ActiveMCS == "El_hin") then
           ! fix MCSmode and change Move.dl, if needed 
         call epdoMixedMCS1
      endif
!        move track tentatively
!              local. coord..

      Move%Track%pos%x = cTrack%pos%x + Move%dl * cTrack%w%x
      Move%Track%pos%y = cTrack%pos%y + Move%dl * cTrack%w%y
      Move%Track%pos%z = cTrack%pos%z + Move%dl * cTrack%w%z

!     *************** see if Move.Track crosses the boundary ******
!              set Cross=t/f, if t, adjust path and set Trunc=t

      call epifCross(lengthToB, icon)
         !  icon =1 , cross. = 0  not cross ; don't care if .not. doMixedMCS
      if( ActiveMCS == "El_hin") then
             ! mixed mode
         if( icon == 0 ) then
           ! not cross the boundary

            call epdoMixedMCS2(icon)
!!!!!!!!!!!!!  remove several lines if no problem
            if(icon == 2 ) then
               write(0,*) ' icon 2 aft mcs2, '
               write(0,*) ' lenB=',lengthToB
               write(0,*) ' MCSmode=', MCSmode
               write(0,*) ' KE,chg=',KEGeV,  cTrack%p%charge
!               stop
            endif
!!!!!!!!!!!!!
            if(icon == 2 ) return  !!!
         endif
      endif

      icon = 0 
!     
!     ************** energy loss consideration ******************
!        if E becomes <= m,  adjust path and set Trunc=t,
!        and reset Cross.
!            kchgPath case Move.dl=0 so light related
!            ptcl will not goto epEloss
      if(cTrack%p%charge .ne. 0 .and. Move%dl .gt. 0.) then
         call epEloss
      else
         Move%dE = 0.
      endif
      if( V1ry == 1 ) then
         Dist2colp = Dist2colp - Move%dl
      endif

!c      endif
      end
!        fix MCSmode and change Move.dl if needed.
      subroutine epdoMixedMCS1
!!!!!!!
      use moddebug
!!!!!!!!!      
      
      use modcMCS
      use modMCScontrol
      implicit none
#include  "Zmedia.h"      
#include "ZepTrackv.h"
      real(8)::u

      if( Move%proc == "hcs" .and. .not. Move%Trunc ) then
               !  ----|---|--
!                    Lh   X
!           case A.  Lh < X 
         pathScm= Move%dl ! (=lHardcm)
         if( KEeV < MCSnow%minNon0mucE ) then
               ! detailed mode, Move.dl unchanged
            MCSmode = "A1"
         else
            MCSmode = "A2"
               !--|---|---|--
!                 t   Lh   X
            call rndc(u)
            tMCScm = lHardcm * u
            Move%dl= tMCScm
         endif
      else
!          case B   X < Lh      !  truncated or not "hcs"
               !--|---|---|--
!                 X   Lh 
         pathScm =Move%dl

         if( KEeV < MCSnow%minNon0mucE ) then ! very unprobable
            MCSmode = "B1"      ! Move.dl unchaged
         else
            call rndc(u)
            tMCScm = lHardcm * u
            if( tMCScm < Move%dl) then  ! X=Move.dl
!            |------|---|---|--
!                t  X   Lh 
               MCSmode = "B2"               
               Move%dl = tMCScm
            else
!            |------|----|---|--
!                   X t  Lh 
               MCSmode = "B3"
            endif
         endif
      endif
      end


      subroutine epdoMixedMCS2(icon)
!!!!!!!!!!!1
      use moddebug
!!!!!!!!!!      
      use modMCScontrol 
      implicit none
#include  "ZmediaLoft.h"      
#include "ZepTrackv.h"
      integer,intent(out):: icon  ! if 2, return immediately

      real(8):: dlNew, A_MCS, muc, uc, cosa, mu
      real(8):: lengthToB
      logical:: dlChange 
      integer:: jcon


!  if hinge is effective, we have to change Move.dl
!  so that it contains total  hinge length
      dlNew = Move%dl
      dlChange=.false.
      icon = 0 

      if( MCSmode == "A1" ) then
         call cfixHardMuc(KEeV, A_MCS, muc, uc)
         uc = 0.      ! for safety for single scattering 
         call csampCSPolA(KEeV, A_MCS, uc, mu, cosa)
         !   embed cosa and adjust p
         call epembedAng(Move%Track, cosa)
         !dl no change
      elseif ( MCSmode == "A2") then
         call ciniFsoft(KEeV, lHardgr) 
         call csampSoftMCSang(mu, cosa)
         !   embed cosa and adjust p
         call epembedAng(Move%Track, cosa)
         ! really move:  tentative Move is copied to cTrack
         call epTraceFE  ! front end for Trace
         cTrack = Move%Track
         ! get distance to the boundary
         call epbndry(Move%boundary, lengthToB, icon)
         if(icon == 2 ) return  !!!!
           ! move streight by lHardcm-tMCScm tentatively
         Move%dl = pathScm - tMCScm

         Move%Track%pos%x = cTrack%pos%x + Move%dl * cTrack%w%x
         Move%Track%pos%y = cTrack%pos%y + Move%dl * cTrack%w%y
         Move%Track%pos%z = cTrack%pos%z + Move%dl * cTrack%w%z
!         see if cross boundary; if cross, adjust Move.Track
         call  epifCross(lengthToB, jcon)
!          independent of jcon
         dlNew = tMCScm + Move%dl
         dlChange = .true.
!          if cross (icon /= 0), proceed to next step
         if( jcon == 0 ) then
           ! if not cross, hard scattering angle should be embeded
            call cfixHardMuc(KEeV, A_MCS, muc, uc)
            call csampCSPolA(KEeV, A_MCS, uc, mu, cosa )
            call epembedAng(Move%Track, cosa)  ! need not copy to Track here
         endif
      elseif( MCSmode == "B2" ) then
         call ciniFsoft(KEeV, lHardgr)
         call csampSoftMCSang(mu, cosa)
                        ! embed cosa
         call epembedAng(Move%Track, cosa)
         ! really move:  tentative Move is copied to cTrack
         call epTraceFE  ! trace info
         cTrack = Move%Track
         ! get distance to the boundary for this direction
         call epbndry(Move%boundary, lengthToB, icon)
         if(icon == 2 ) return  !!!!
           ! move streight by X-tMCScm tentatively

         Move%dl = pathScm - tMCScm

         Move%Track%pos%x = cTrack%pos%x + Move%dl * cTrack%w%x
         Move%Track%pos%y = cTrack%pos%y + Move%dl * cTrack%w%y
         Move%Track%pos%z = cTrack%pos%z + Move%dl * cTrack%w%z
!         see if cross boundar; if cross, adjust Move.Track
         call  epifCross(lengthToB, jcon)
!          if not cross (icon == 0), proceed to next step
!           nothing to do here. Move.dl  has been modifed if cross.
         dlNew = tMCScm + Move%dl
         dlChange = .true.
      elseif( MCSmode == "B1" ) then
               ! nothing to do
      elseif( MCSmode == "B3" ) then
               ! nothing to do
      else
         write(0,*) ' MCSmode =', MCSmode, 'invalid'
         stop
      endif
      if( dlChange) then
         Move%dl = dlNew
         Move%dx = Move%dl/Media(MediaNo)%gtocm*Media(MediaNo)%rhoc
         Move%dt = Move%dx*Media(MediaNo)%X0g ! in r.l
      endif
      end

      subroutine epembedAng(aTrack, cosa)
      implicit none
#include  "Zmedia.h"      
#include "ZepTrackv.h"
      
      real(8):: tmp, sina, cs, sn
       type(epTrack)::  aTrack
      real(8),intent(in):: cosa

       type(epDirec)::  w

      tmp=max(1.d0-cosa*cosa, 0.d0)
      sina = sqrt(tmp)
        ! azimuthal ang.
      call kcossn(cs,sn)
      w%x = cs*sina
      w%y = sn*sina
      w%z = cosa
      call eptransVect(aTrack%w,  w, aTrack%w)
!        energy unchaged;   modify p
      call epe2p(aTrack)
      end

      subroutine epqProc(proc)
!          inquire the current process fixed
      implicit none
#include  "ZepManager.h"
#include  "Zmedia.h"
#include  "ZepTrackv.h"
      character(8),intent(out)::proc
      proc = Move%proc
      end
!     *******************************
      subroutine epCross(icon)
      implicit none
!     manager when a particle crosses the boudarynon
#include "ZmediaLoft.h"      
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
#include "Zcode.h"

      integer icon ! input/ouptut  1--> discard this particle
!                                 other values are not given here
      integer cnx, info
       type(epPos):: postemp
       type(epDirec)::  dirtemp

       type(epTrack)::  aTrack
      integer n

!
!      save the cTrack
      aTrack = cTrack
!       cross the boundary;  use  new angle at the boundary.
!      (not cTrack.w.x etc)

      cTrack%pos%x = Move%boundary%x +
     *     EpsLeng* Move%Track%w%x  
!cccc          if EpsLeng and w.x are both small, and x is large
!cccc          no move may happen, avoid such case.
!cccc               You may give larger EpsLeng in input data.
!cccc
!cc      if( cTrack.w.x .ne. 0. ) then
!cc         n = 2
!cc         do while( cTrack.pos.x .eq. Move.boundary.x .and.
!cc     *             n .lt. 10 ) 
!cc           cTrack.pos.x = Move.boundary.x +
!cc     *           n*EpsLeng* Move.Track.w.x
!cc            n = n * 2
!cc         enddo
!cc      endif
!cccccc
      cTrack%pos%y = Move%boundary%y +
     *     EpsLeng* Move%Track%w%y
!cc      if( cTrack.w.y .ne. 0.) then
!cc         n = 2
!cc         do while( cTrack.pos.y .eq. Move.boundary.y
!cc     *             .and. n .lt. 10 ) 
!cc           cTrack.pos.y = Move.boundary.y +
!cc     *           n*EpsLeng*  Move.Track.w.y
!cc           n = n * 2
!cc         enddo
!cc      endif
!cc
      cTrack%pos%z = Move%boundary%z +
     *     EpsLeng* Move%Track%w%z
!cc       if(cTrack.w.z .ne. 0. ) then
!cc         n = 2
!cc         do while( cTrack.pos.z .eq. Move.boundary.z  .and.
!cc      *          n .lt. 10 )
!cc            cTrack.pos.z = Move.boundary.z +
!cc      *           n*EpsLeng* Move.Track.w.z
!cc            n = n * 2
!cc          enddo
!cc       endif
!cc 
#ifdef  SUBSTREC
      cTrack%w = Move%Track%w
#else
      call epsubvec(Move%Track%w, cTrack%w)
#endif

!         new comp. number;  new position would be, in principle,
!         in a new component. but scattering at the boundary may
!         cause the new pos to be in the same component.
      call eppos2cn(Cn, cTrack, cnx)   ! cnx is not yet set in cTrack
 

      if(Cn .eq. cnx) then
         Move%Cross = .false.
      elseif(cnx .gt. Det%nct) then
         info = 0   ! exiting to void
      else
         info = -cnx   ! exiting to cnx
      endif
      cTrack%cn = cnx
!         update coord. to local one in new  comp.; this part was inside
!              " elseif(Move.Cross) then"
!         which appears later; since it is used in LightAbBndry
      if(Move%Cross) then
         call epl2w(Cn, cTrack%pos, postemp)
         call epw2l(cnx, postemp,  cTrack%pos)
!!!   call epl2wd(Cn, cTrack%w,  dirtemp)   ! 2016 Sep
          !  replaced by next
         call epl2wdm(Cn, cTrack%w,  dirtemp, cTrack%p)   ! 2016 Sep
!cc         call epw2ld(cnx, dirtemp, cTrack.w)
         call epw2ldm(cnx, dirtemp, cTrack%w, cTrack%p)
      endif

       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>light
      if(Move%Cross .and.  cTrack%p%code < 0 ) then
         call epLightAtBndry( cnx, icon )
            ! cTrack, Move and Cn are implicit  in/out param.
            ! if refraction happens, cTrack's direction becomes 
            !    refracted direc. pos, cnx are unchaged
            !    icon = 0
            ! if reflection happens, Move.Cross=F, cnx will become Cn.
            !    cTrack's pos  is made to be that of  Move. 
            !    icon = 0
            ! if absorbed at the boundary,  icon = 1 Move.Cross=F
            ! if light passes thru the component, icon = 0. cTrack
            !   unchaged. cnx unchaged
      endif
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<
      if(Det%cmp(Cn)%CountIO .ge. 2 .and. Move%Cross) then
!          user hook for  counting exiting ptcls 
!                    exiting from Cn to cnx
         if(Light > 0 .and. aTrack%p%code == klight) then
            call epLightPC(info)  ! photon counter  for exiting light
         endif
         call userbd(info, aTrack, Move, Media(MediaNo))
         if(Move%Abort .ne. 0) then
            if(Move%Abort .eq. 3) then
               icon = 1         ! discard this particle
               Move%Abort = 0   ! but continue simulation
            else
               call  epempty    ! empty the sack. discard ev.
               call  epSkipUpdateNo
            endif
         endif
      endif


      if(cnx .gt. Det%nct) then
!          void
         icon = 1
      elseif(Move%Cross) then
         info = Cn            ! save the current cn
!            update comp. info (Cn, MediaNo etc are updated)
         call epnewComp(cTrack)
#ifdef  SUBSTREC
         aTrack = Move%Track
#else
         call epsubstTRK(aTrack, Move%Track)
#endif
#ifdef SUBSTREC         
         Move%Track = cTrack      ! copy the track which is  now in new comp.
#else
         call epsubstTRK(Move%Track, cTrack)
#endif

         if(Det%cmp(Cn)%CountIO .eq. 1 .or.
     *      Det%cmp(Cn)%CountIO .eq. 3 ) then
!             count entering ptcls;  from info to Move.Track.cn (=Cn)
            if( Det%cmp(Cn)%CountIO == 1 ) then
               if(Light > 0  .and. aTrack%p%code == klight) then
                  call epLightPC(info) ! photon counter  for entering light
               endif
            endif
            call userbd(info, aTrack, Move, Media(MediaNo)) 
            if(Move%Abort .ne. 0) then
               if(Move%Abort .eq. 3) then
                  icon = 1      ! discard this particle
                  Move%Abort = 0 ! but continue simulation
               else
                  call  epempty ! empty the sack. discard ev.
                  call  epSkipUpdateNo
               endif
            endif
         endif
      endif                  
      end
!         for ibm
      subroutine epsubvec(inp,out)
      implicit none
#include "Zep3Vec.h"
       type(ep3Vec)::   inp,out
      out = inp
      end
      subroutine epsubstTRK(left, right)
      implicit none
#include "ZepTrack.h"
       type(epTrack):: left, right
      left = right
      end

!     **********************
      subroutine epNonEleMag(media)
!               obso      
!         media in modXsecMedia is now renamed as xmedia
!       element    //                             xelement  
!       NoOfMedia  //                          as notused
!      use modXsecMedia, xmedia=>media, xelement=>element,
!     *  dontuse=>NoOfMedia
!               new 

      implicit none
#include "Zmedia.h"
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcode.h"

      type(epmedia),intent(in):: media
!////////////////
!      real(8):: pathp, pathb, pathn  
!///////////       
      integer jcon
!          sample interaction length for non-e/g
      real*8 prob, path, tokgpm2

!          set cosmos condition(TrackBefMove,
      tokgpm2 = media%X0g *10.d0 ! r.l --> kg/m^2

!      call ep2cosCond
      call cfixModel( cTrack%p )
!!      call ciniSmpIntL --->
!!      call ciniIntInf  !         moved to before process sampling
      call  cepSampNEPIntL(media, cTrack%p)
      
!      call epsmpNEPIntL( media )
!         reset muon interaction conditions

      call epfixProc(3, media%rho*media%rhoc, path) ! path is in  m

      if( Move%proc == 'coll' ) then
         call  cseeColPossible( cTrack%p, jcon)
         if(jcon == -1) then
            Move%proc = 'decay'
            call cresetIntInf   ! this rest is needed
         endif
      endif

      Move%dl = path*1.d2       ! cm
      Move%dx = Move%dl * media%rhoc / media%gtocm
      Move%dt = Move%dx/media%X0g ! in r.l
!      Move%dl = Move%dx * media%gtocm /  media%rhoc  ! in cm
      end
!     ******************
      subroutine epEloss
      use modMCScontrol
!!!!!!!!!
      use moddebug
!!!!!!!!!      
      use moddedx
      use modGUI
      implicit none
#include "ZepTrackp.h"
#include "ZmediaLoft.h"      
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "Zcode.h"

      type(gui),save:: sinfo  ! for quenching
      integer::modif
      integer k
      
      real*8  dx, up, dedxmu, dedlmu,
     *      s1, s2, scol2, schg2, sigma
      real*8 cupsilon, csyncTELoss, cf, qf, cff
      real(8):: dEdxbyBrems, dEdxbyBremsf
      real(8):: dEdlBrems, FbyR
      logical lowehvyion
      real(8),save::  dedl
      real(8):: dedxout, dedlout
       !>>>>>>>>>>>>>>>light
      integer::epLightGetd  !  function for extracting d part from countDE
       !<<<<<<<<<<<<<<<<<
      qf = 0.  ! quenching fraction. 
      k = cTrack%p%code
!              compute energy loss rate
      if(EdepdEdx) then
         dEdxbyBrems = 0.
         dEdxbyBremsf = 0.
         dEdlBrems = 0.
         lowehvyion  = .false.
         if(cTrack%p%fm%p(4) .le. cTrack%p%mass) then
            dedx = 0.
            Move%dE = 0.
            Move%dEeff = 0.
            Move%dEioni = 0.
         else
            if(k .eq. kelec) then
!                            dedx is in GeV/(g/cm^2)
               call epdedxe(Media(MediaNo), cTrack%p, dedx, dedxf)
               dedxmu = 0.
            else
               if( abs( cTrack%p%charge ) .gt. 1 ) then
                  if( k /= kgnuc ) then  
                       ! may appear once 1/10^8 col. ?
                     call epStrange
                     dedx = 0.
                     Move%dE = 0.
                     Move%dEeff = 0.
                     Move%dEioni = 0.
                     return
                  endif
                  lowehvyion  = 
     *              (cTrack%p%fm%p(4)-cTrack%p%mass)/cTrack%p%subcode
     *               < 0.7
                  call epdedxhvy(Media(MediaNo), cTrack%p,
     *               dedx,dedxf)
               else
                  call epdedxNone(Media(MediaNo), cTrack%p, dedx,dedxf)
                  call epdedxTargetBrems(Media(MediaNo), cTrack%p, 
     *               dEdxbyBrems, dEdxbyBremsf)
               endif

               if(cTrack%p%code .eq. kmuon) then
                  call epmudEdx(Media(MediaNo)%mu%MuNI, 
     *                 Media(MediaNo)%mu%MuBr, 
     *                 Media(MediaNo)%mu%MuPr, 
     *                 Media(MediaNo),
     *                 cTrack%p%fm%p(4), dedxmu)
               else
                  dedxmu = 0.
               endif
            endif

            if( abs(cTrack%p%charge) .gt. 1 ) then
!                 below for test; use above  normally
!            if( cTrack.p.charge >= 1) then ! ????
               modif = Det%cmp(Cn)%modifier
               sinfo%modif = modif
!                   next Birks is for Talre, Birks or Log
               if(modif > 0 .or.
     *            Media(MediaNo)%Birks /= ' ') then
!                even if no quenching is specified in media file,
!                modifier can  force quenching. inside of the next,
!                kind is checked if really quenching is specified by modifier

                  if( HowQuench == 0 )  then
!                   employ cf as it is (same as old) "dedx" 
                     call epOrgCorrec( modif,
     *                    Media(MediaNo), cTrack%p, dedx, cf)
                     qf= 1.0
                  elseif( HowQuench >= 1 ) then
                     call epOrgCorrec( modif,
     *                    Media(MediaNo), cTrack%p, dedxf, cff)
                     sinfo%dedx = dedx
                     sinfo%dedxf = dedxf
                     call epGUI(1, sinfo)
                     cf = sinfo%cf
                     qf = sinfo%qf
                  else
                     write(0,*) ' HowQuench=',HowQuench,
     *               ' invalid '
                     stop
                  endif
               else
                  cf = 1.
               endif
            else
               cf = 1.
            endif

            dedl = dedx /Media(MediaNo)%gtocm * ! GeV/cm
     *         Media(MediaNo)%rhoc
            dedlmu = dedxmu /Media(MediaNo)%gtocm * ! GeV/cm
     *         Media(MediaNo)%rhoc
            if( btest(TargetElecBrems, 2)) then
!                         this will overestimate for thin media
               dEdlBrems = dEdxbyBremsf/Media(MediaNo)%gtocm * 
     *         Media(MediaNo)%rhoc   
            elseif( btest(TargetElecBrems, 0) .or.
     *           btest(TargetElecBrems, 1)  ) then
!                     this will underestimate
               dEdlBrems = dEdxbyBrems/Media(MediaNo)%gtocm * 
     *         Media(MediaNo)%rhoc   
            endif
            if( abs( epLightGetd( Det%cmp(Cn)%CountDE) ) .eq. 2) then
!                Enery loss fluctuation; use  Urban model for singl charge
!                 or  high energy
               if( lowehvyion ) then
                  call epdedxflhv(Media(MediaNo), cTrack%p,
     *                scol2, schg2)
                  sigma = sqrt( (scol2 + schg2)* Move%dx )
                  call kgauss(0.d0, sigma, s1, s2)
                  Move%dEioni = max(dedx*Move%dx + s1, 0.d0)
               else
                  call epUrban(Media(MediaNo)%urb, dedl, Move%dl,
     *                 cTrack%p,  Move%dEioni)
               endif
!                    Move.dEioni = GeV  in 'dl' cm
            else
               Move%dEioni = dedl * Move%dl
            endif
            Move%dE = Move%dEioni + dedlmu*Move%dl  
     *            + dEdlBrems*Move%dl  ! total loss without quenching
            Move%dEeff = Move%dEioni *(qf* cf + (1-qf))
     *            + dedlmu*Move%dl  
     *            + dEdlBrems*Move%dl   ! effective loss with quenching
         endif
      else
!            use constant energy loss at minimum ionization; actually not used
!            In this case, knock-on process should be prohibitted  by
!            giving large Tcut.   Birks etc correction not applied.
         dedx = Media(MediaNo)%dEdxatp3m * cTrack%p%charge**2  ! GeV/(g/cm2)
         Move%dE = dedx * Move%dx
         Move%dEioni = Move%dE
         Move%dEeff = Move%dEioni
         SumDe = SumDe + Move%dE 
      endif
!    ****************** synchrotron loss; only in 'sp' and for electrons
!        This is normally completely negligible
      if(Sync .eq. 1) then
         if(Media(MediaNo)%name .eq. 'sp' ) then
            if(MagField .gt. 0 .and. k .eq. kelec) then
               up = cupsilon(cTrack%p, Bfield) ! Upsilon value
               dedl = csyncTELoss(up) * 10.d-2
                               !  GeV/cm: dedl by cosmos is GeV/m
               Move%dE = Move%dE + dedl*Move%dl
            endif
         endif
      endif
!
!  
      Move%Track%p%fm%p(4) = cTrack%p%fm%p(4) - Move%dE
      if(Move%Track%p%fm%p(4) .le. Move%Track%p%mass) then
!          this will not happen if muon specific loss exists.
         Move%dE = max(cTrack%p%fm%p(4)- cTrack%p%mass, 0.d0)

         if(dedx .le. 0. .or. Move%dE .eq. 0.) then
            dx = 0.
            Move%Cross = .false.
            Move%Trunc = .false.
            Move%dE = 0.
            Move%dEeff = 0.
            Move%dEioni = 0.
         else
            dx = max(Move%dE/dedx, 0.d0)
            Move%Trunc=.true.
            Move%Cross = .false.   
            Move%dEeff = Move%dE*(cf*qf+(1.0-qf))
            Move%dEioni = Move%dE
         endif

         Move%Track%p%fm%p(4) = Move%Track%p%mass
         Move%Track%p%fm%p(1:3) = 0.
!           in priciple, min need not be take
!            but is some strange case, dl becomes large
!            than old dl. 
         Move%dl =min( dx * Media(MediaNo)%gtocm /
     *     Media(MediaNo)%rhoc, Move%dl)
         Move%dx = dx
!!!!!!!!!!Jan. 2017.
! in the case of El_hin, Move%dl is long (
!     from previous hard cs point, while cTrack point is
!     hinge point, so Moved point below could be outside
!     of the prsent volume. Althouth the particle is already
!     dead, in the next step, it is treated as if a live
!     particle and the boundary is searched for. If the
!     position is not in the true component, boundary seachh
!     error happens. It can be recovered automatically,but
!     if it happens > 10 times /event, stop is forced.
!     To avoid such, we keep the moved particle position
!     at the current hinge point. The energy loss is
!     in the correct compoenent, so little difference
!     of stopping point dosn't matter.
         if(ActiveMCS == "El_hin") then
            Move%dl=0.
            MOve%dx =0.
         endif
!!!!!!!!!
         Move%Track%pos%x = cTrack%pos%x + Move%dl* cTrack%w%x 
         Move%Track%pos%y = cTrack%pos%y + Move%dl* cTrack%w%y
         Move%Track%pos%z = cTrack%pos%z + Move%dl* cTrack%w%z
      endif
      SumDe = SumDe + Move%dE
      end
!     ****************************************
      subroutine epprog(media)
      use modEMcontrol
      implicit none
#include  "Zmedia.h"    
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zevhnv.h"
      type(epmedia),intent(in):: media

      rrho = media%rhoc
      
      call cepSampGIntL( media, cTrack%p)
      if( MagPair == 1 ) then
         call eppos2B(cTrack, Bfield)
         call epmpairp(cTrack%p, Bfield, Xai, pairmfp, dl) ! dl in m
! for comparison with other processes, convert it to r.l
         ! Xai is saved in modEMcontrol
         dl = dl * 100.d0/ media%X0 
         call csetIntInf(dl, .false., 'mpair')
      endif

      call epfixProc(1,  media, path)
      Move%dl = path * 1.0d2
!            Move%dl = Move%dx * meddia%gtocm / media%rhoc, so
      Move%dx = Move%dl * media%rhoc / media%gtocm
!        Move%dx = Move%dt * mediaX0g, so  
      Move%dt = Move%dx/media%X0 ! in r.l


      end
!!!      subroutine epprog
!!!       use modXsecMedia, xmedia=>media, xelement=>element,
!!!     *     dontuse=>NoOfMedia
!!!       implicit none
!!!!
!!!!     5     processes, i.e, pair creation, compton, Photo-electric
!!!!     effect, coherent scatt.
!!!!      and Photo-production of pions are considered.
!!!!
!!!!
!!!#include  "ZepTrackv.h"
!!!#include  "ZepTrackp.h"
!!!#include  "Zevhnv.h"
!!!
!!!      real*8 tcomp, tphot, tpair, tgp, tcoh,  t
!!!      real*8 E, prob, xprob(5), txray(5)
!!!      real*8 xs, mfp
!!!      real*8  pairmfp, dl, dx, tmpair, u
!!!      real*4  xsec(5)  !  coh, incoh,  P.E  1/(g/cm2)
!!!!      real*8  Excom1, Excom2  !now in ZepTrackp.h
!!!      integer icon
!!!!            where xcom data is used. 
!!!!       Excom1: compton/photo abs/coherent scat
!!!!       Excom2: pair; default is  not use xcom
!!!!              both must be < 100 GeV
!!!!      data Excom1/1.d-3/, Excom2/1.d-3/
!!!!      save Excom1, Excom2
!!!
!!!      E = cTrack%p%fm%p(4)
!!!      if(E .le. EupperBndCS) then
!!!         xprob(3)= 0.
!!!         if( Media(MediaNo)%xcom%size .gt. 0 .and.
!!!     *        E .lt. Excom1 ) then
!!!!               below 1MeV, use accurate xs.
!!!            call epXrayp(Media(MediaNo), E, 1,  3,  xprob, txray)
!!!            tcomp=txray(2)
!!!            tcoh = txray(1)
!!!            tphot = txray(3)
!!!         else
!!!            tcoh = 1.e35
!!!            call epcompp(Media(MediaNo), E, prob, tcomp)
!!!         endif
!!!         if(Photo) then
!!!            if(xprob(3)  .eq. 0.) then
!!!               call epphotoEp(Media(MediaNo),  E, prob, tphot) ! v8.0
!!!            endif
!!!         else
!!!            tphot = 1.e35
!!!         endif
!!!      else
!!!         tcomp =1.e35
!!!         tphot= 1.e35
!!!         tcoh = 1.e35
!!!      endif
!!!      if(E .gt. ElowerBndPair) then
!!!         if( Media(MediaNo)%xcom%size .gt. 0 .and.
!!!     *        E .lt. Excom2 ) then
!!!            call epXrayp(Media(MediaNo), E, 4, 5,  xprob, txray)
!!!            prob= xprob(4)+xprob(5)
!!!            call rndc(u)
!!!            tpair = -log(u)/prob
!!!         else
!!!            call epPrSampP(Media(MediaNo), E, prob, tpair)
!!!         endif
!!!      else
!!!         tpair=1.e35
!!!      endif
!!!      if(IncGp > 0 .and. E .gt. 153.d-3) then ! > 0.152431798977028
!!!         call ep2cosCond
!!!         call cfixModel( cTrack%p )
!!!         call cgetPhotoPxs(ActiveMdl2, cTrack%p, xmedia(mediumNo),
!!!     *      xs, mfp)
!!!         prob = xs*Media(MediaNo)%mbtoPX0 ! prob/r.l
!!!         call rndc(u)
!!!         tgp = -log(u)/prob          ! sampled path in r.l
!!!      else
!!!         tgp=1.e35
!!!      endif
!!!
!!!      if(MagPair .eq. 1) then
!!!         call epmpairp(cTrack%p, Bfield, Xai, pairmfp, dl)
!!!         dx = dl / Media(MediaNo)%gtocm *
!!!     *     Media(MediaNo)%rhoc
!!!         tmpair = dx / Media(MediaNo)%X0g
!!!      else
!!!         tmpair = 1.e35
!!!      endif
!!!      t=min(tpair, tcomp, tphot, tgp, tmpair, tcoh)
!!!      if(t .eq. tpair) then
!!!         Move%proc='pair'
!!!      elseif(t .eq. tcomp) then
!!!         Move%proc='comp'
!!!      elseif(t .eq. tphot) then
!!!         Move%proc='phot'
!!!      elseif(t .eq. tcoh) then
!!!         Move%proc='coh'
!!!      elseif(t .eq. tgp) then
!!!         Move%proc='photop'
!!!      else
!!!         Move%proc='mpair'
!!!      endif
!!!      Move%dt = t   ! in r.l
!!!      Move%dx = Move%dt * Media(MediaNo)%X0g
!!!      Move%dl = Move%dx * Media(MediaNo)%gtocm /
!!!     *     Media(MediaNo)%rhoc
!!!   end
      subroutine  epproe(media)
      use modEMcontrol
      implicit none
#include  "Zmedia.h"
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"

      type(epmedia),intent(in):: media

      rrho = media%rhoc
      
      call cepSampEIntL(media,   cTrack%p)
      call epMCShard      
!!!!                     only if X0 > 10 km; may be 30 km is o.k
      if(Sync .eq.  2 .and.
     *   media%X0/media%rhoc .gt. 10.d5 ) then
!     sample synchrotron emission path
         call eppos2B(cTrack, Bfield)
         call epsyncp(cTrack%p, Bfield, Upsilon, syncmfp, dl) !dl in m
          ! Upsilon is saved in modEMcontrol
         dx = dl / media%gtocm *  media%rhoc
         dt = dx / media%X0g
         call csetIntInf(dt, .false., 'sync')
      endif

      call epfixProc(2,  media, path) ! path is in m

      Move%dl = path * 1.0d2    ! cm
!            Move%dl = Move%dx * meddia%gtocm / media%rhoc, so
      Move%dx = Move%dl * media%rhoc / media%gtocm   ! g/cm2
!        Move%dx = Move%dt * mediaX0g, so  
      Move%dt = Move%dx/media%X0 ! in r.l

            
      end      
!!!!     *************
!!!      subroutine  epproe
!!!!     *************
!!!!!!!!!!!!
!!!      use moddebug
!!!!!!!!!!!!      
!!!      use modMCScontrol
!!!      implicit none
!!!!
!!!!   electron:     brems, knock-on, and anihilation, synchrotron
!!!!                 radiation are considered
!!!!
!!!#include  "ZepTrackv.h"
!!!#include  "ZepTrackp.h"
!!!!
!!!      real*8 E, prob,  tbrem, tknock, tanihi, t, dt, dl, dx,
!!!     *   syncmfp
!!!      real(8)::u
!!!      logical:: posok, numok
!!!      type(epPos):: posw
!!!      integer:: i
!!!       
!!!      E = cTrack%p%fm%p(4)
!!!!             sample path for brems
!!!      call epBrSampP(Media(MediaNo),  E, prob, tbrem)
!!!
!!!      if(Knckon) then
!!!         if(cTrack%p%charge .eq. -1) then
!!!            call epmollerp(Media(MediaNo), E, RecoilKEmin, prob, tknock)
!!!         else
!!!            call epbhabhap(Media(MediaNo), E, RecoilKEmin, prob, tknock)
!!!         endif
!!!         if(tbrem .le. tknock) then
!!!            t = tbrem
!!!            Move%proc='brem'
!!!         else
!!!            t = tknock
!!!            Move%proc='knoc'
!!!         endif
!!!      else
!!!         t = tbrem
!!!         Move%proc='brem'
!!!      endif
!!!      if(cTrack%p%charge .eq. 1 .and. E .lt. Eanihi) then
!!!         call epanihip(Media(MediaNo), E, prob, tanihi)
!!!         if(tanihi .lt. t) then
!!!            t = tanihi
!!!            Move%proc='anih'
!!!         endif
!!!      endif
!!!
!!!      if( doNewMCS ) then
!!!         if( (Media(MediaNo)%name /= 'sp' ) .and.
!!!     *       (Media(MediaNo)%name /= 'hollow')  ) then
!!!            call cfixMCSmodel( cTrack%p ) ! see energetically OK ?
!!!
!!!            if( ActiveMCS /= 'Mol') then
!!!
!!!               if( MCSzCond > 0 ) then
!!!!                see if z is in the specifed region; else -> Mol
!!!                  posok = .false.
!!!                  call epl2w(cTrack%cn, cTrack%pos, posw)
!!!                  do i = 1, MCSzCond
!!!                     if(real(MCSzRange(i)) <= posw%z ) then
!!!                        if( posw%z <= imag(MCSzRange(i)) ) then
!!!                           posok = .true.
!!!                           exit
!!!                        endif
!!!                     endif
!!!                  enddo
!!!               else
!!!                  ! no check about Z
!!!                  posok=.true.
!!!               endif
!!!
!!!               if( MCSnumCond > 0 ) then
!!!                  if( MCSdebug ) then
!!!                     call epl2w(cTrack%cn, cTrack%pos, posw)
!!!                  endif
!!!                  ! see if cn is in the ranage
!!!                  numok = .false.
!!!                  do i = 1, MCSnumCond
!!!                     if(MCSnumRangeMin(i) <=  cTrack%cn ) then
!!!                        if( cTrack%cn <= MCSnumRangeMax(i) ) then
!!!                           numok = .true.
!!!                           exit
!!!                        endif
!!!                     endif
!!!                  enddo
!!!               else
!!!                  ! no check about cn
!!!                  numok = .true.
!!!               endif
!!!
!!!
!!!               if( MCSandor == 'and' ) then
!!!                  if( posok .and. numok ) then
!!!                     if( MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!                        !  keep current one
!!!                     endif
!!!                  else
!!!                     if( .not. MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!!                      keep current
!!!                     endif
!!!                  endif
!!!                  
!!!               else             ! 'or'
!!!                  if( posok .or.  numok ) then
!!!                     if(  MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!!                         keep current
!!!                     endif
!!!                  else
!!!                     if( .not. MCSrevert ) then
!!!                        ActiveMCS = 'Mol'
!!!                     else
!!!!                          keep current
!!!                     endif
!!!                  endif
!!!               endif
!!!               if( MCSdebug ) then 
!!!                  if( ActiveMCS == 'Mol') then
!!!                     write(0,'(a, 1p,3g14.4)') 'mol: ',  posw
!!!                  else
!!!                     write(0,'(a, 1p,3g14.4)') 'hin: ',  posw
!!!                  endif
!!!               endif
!!!            else
!!!               ActiveMCS = 'Mol'
!!!            endif
!!!         else
!!!            ActiveMCS = 'Mol'
!!!         endif
!!!      endif
!!!      
!!!      if( ActiveMCS /= 'Mol' ) then
!!!         call cfixMixedConst(MediaNo, int(cTrack%p%charge))
!!!      endif
!!!
!!!
!!!      if( ActiveMCS == 'El_hin') then
!!!             ! sample hard scattering  mfp g/cm2.
!!!         call cgetLamh( KEeV, mfpHardgr)  ! mfp of hard cs
!!!!  in g/cm2
!!!         call rndc(u)
!!!         lHardgr = -log(u)* mfpHardgr  ! g/cm2
!!!         lHardrl =lHardgr / Media(MediaNo)%X0g    ! rl
!!!          !  r.l
!!!          ! next is  cm length. this might  be used later for soft
!!!          ! mcs treatment.
!!!         lHardcm = lHardgr* Media(MediaNo)%gtocm
!!!     *        / Media(MediaNo)%rhoc
!!!         if( lHardrl < t ) then
!!!            Move%proc="hcs"  ! hard Coulomb scattring
!!!            t = lHardrl
!!!         endif
!!!      endif
!!!
!!!      Move%dt = t
!!!      Move%dx = Move%dt * Media(MediaNo)%X0g
!!!      Move%dl = Move%dx * Media(MediaNo)%gtocm /
!!!     *    Media(MediaNo)%rhoc
!!!
!!!!                     only if X0 > 10 km; may be 30 km is o.k
!!!      if(Sync .eq.  2 .and.
!!!     *    Media(MediaNo)%X0/Media(MediaNo)%rhoc .gt. 10.d5 ) then
!!!!          sample synchrotron emission path
!!!         call epsyncp(cTrack%p, Bfield, Upsilon, syncmfp, dl)
!!!         dx = dl / Media(MediaNo)%gtocm *
!!!     *    Media(MediaNo)%rhoc
!!!         dt = dx / Media(MediaNo)%X0g
!!!         if(dt .lt.  t) then
!!!            Move%dt = dt
!!!            Move%dx = dx
!!!            Move%dl = dl
!!!            Move%proc = 'sync'
!!!         endif
!!!!      else
!!!      endif
!!!      end
!     ******************
      subroutine eptrunc
      use modV1ry
      implicit none
!          truncate path if it is too long.  
!          Move.dt,dx,dl is adjusted if truncated, and Trunc is set to T.
!
#include  "ZmediaLoft.h"          
#include  "ZepTrackv.h"
#include  "ZepTrackp.h"
#include  "Zcode.h"

!           max length movable
      real*8 tmax, r
      if(cTrack%p%charge .eq. 0) then
!           neutral particle can move any length basically
         Move%Trunc = .false.
!           but we set max r.l be 100 r.l so that
!           boundary calculation has high accuracy
         if(Move%dt .gt. 100.) then
            Move%dt = 100.
            Move%dx = Move%dt * Media(MediaNo)%X0g
            Move%dl = Move%dx * Media(MediaNo)%gtocm /
     *           Media(MediaNo)%rhoc
            Move%Trunc = .true.
         endif            
      elseif( Light > 0  .and.  cTrack%p%code == kchgPath ) then
         Move%Trunc = .false.   
             ! Move.dt  etc is 0

!<<<<<<<<<<<<<<<<<
!      elseif(Move.proc .eq. 'no') then
!         Move.Trunc=.false.
!         if(MagField .gt. 0 .or.  ElecField .gt. 0) then
!            Move.dt = min( max( Tcoefx*cTrack.p.fm.p(4), Tminx), 
!     *           Move.dt, MaxPath )
!         endif
      else

         tmax = min(
     *     max( Tcoefx*(cTrack%p%fm%p(4)-cTrack%p%mass), Tminx), 
     *     MaxPath )



         if(Move%dt .lt. tmax) then
            Move%Trunc=.false.
         else
            Move%Trunc = .true.
            Move%dt = tmax
            Move%dx = Move%dt * Media(MediaNo)%X0g
            Move%dl = Move%dx * Media(MediaNo)%gtocm /
     *           Media(MediaNo)%rhoc
         endif
         if(MagField .gt. 0 ) then
            call epGetB
            call epmagDefR(cTrack, Bfield, r)
!              assume we can go r/10 streight way.
            if(Move%dl .gt. r/10.d0) then
               Move%dl = r/10.d0
               Move%dx = Move%dl/Media(MediaNo)%gtocm *
     *          Media(MediaNo)%rhoc
               Move%dt = Move%dx/ Media(MediaNo)%X0g
               Move%Trunc = .true.
            endif
         endif
         if(ElecField .gt. 0) then
            call epGetE
!                 to be implemented in future
         endif  
      endif
      if( V1ry == 1 .and. .not. Move%Trunc) then
         call eptruncV1ry
      endif
      
      end
!     ****************
      subroutine epGetB
      implicit none
#include  "ZmediaLoft.h"          
#include "ZepTrackv.h"
#include "ZepTrackp.h"
!          if MagField=1, Bfield is unchanged from the start time
      if(MagField .eq. 2) then
         call eppos2B(cTrack, Bfield)
!             Bfield is in local coordinate
      endif
      end
!     *****************
      subroutine epGetE
      implicit none
#include  "Zmedia.h"          
#include "ZepTrackv.h"      
#include "ZepTrackp.h"      
!          if ElecField=1, Efield is unchanged from the start time
      if(ElecField .eq. 2) then
         call eppos2E(cTrack, Efield)
      endif
      end
!     *******************
      subroutine epifCross(el, kcon )
      implicit none
#include "ZmediaLoft.h"
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
!             see if Move.Track crosses a boundary and set Cross = t/f
!        If t, reset Move.Track.pos just before the boundary.
!         kcon  = 0--> not cross
!               = 1--> cross
      integer kcon, icon


      real*8 el

      if(Move%dl .lt. el) then
         kcon = 0
         Move%Cross = .false.
      else
         kcon = 1
         Move%Cross = .true.
         Move%Trunc = .true.

         Move%Track%pos%x = Move%boundary%x -
     *        EpsLeng* cTrack%w%x
         Move%Track%pos%y = Move%boundary%y -
     *        EpsLeng* cTrack%w%y
         Move%Track%pos%z = Move%boundary%z -
     *        EpsLeng* cTrack%w%z
         Move%dl = el - EpsLeng
!                X0=X0g/rho
         Move%dt = Move%dl/Media(MediaNo)%X0 *
     *         Media(MediaNo)%rhoc
!            gtocm = X0/X0g
         Move%dx = Move%dl/Media(MediaNo)%gtocm *
     *         Media(MediaNo)%rhoc
      endif
      end
!
      subroutine epchckE0(aTrack, icon)
      use epModify
      use moddedx
      implicit none
#include "ZmediaLoft.h"          
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
#include "Zcode.h"
#include "Zmass.h"
!            aTrack is examined if it's  energy is too low
!
       type(epTrack)::  aTrack   ! input. this track's energy is examined
      integer icon              ! output. 0. the particle is still alive
                                !         1. death. 

      logical ok, kbtest, needdedx
      integer k
      real*8 ke
      integer,save::nstrange=0
!////////////////                                                              
!      logical show
!      common /showshow/ show
!///////////       

      ke = aTrack%p%fm%p(4)- aTrack%p%mass
      k =  aTrack%p%code
! >>>>>>>>>>>>>>>>>>light
      if(k .eq. klight) then
         if(Light <=  0 ) then
            write(0,*) ' Light =0 but light appeared'
            icon = 1   ! light is not treated, so death
            stop
         endif
         call epLightchkE( aTrack, icon) ! check wave length
         return  ! *************
      endif
!<<<<<<<<<<<<<<<<<<<<light
      if(k .eq. kelec) then
          if(aTrack%p%charge .eq. -1) then
             ok=aTrack%p%fm%p(4) .gt. EminElec
          else
!              positron
             ok=aTrack%p%fm%p(4) .gt. EminGamma 
!             in very rare case, ke becomes ~ -10^-6
!             so force to 0. v9.17
             if( ke < 0. ) then
                aTrack%p%fm%p(4)=aTrack%p%mass
                aTrack%p%fm%p(1:3) = 0.
             endif   
          endif
       elseif(k .eq. kphoton) then
          ok= aTrack%p%fm%p(4) .gt. EminGamma
       elseif(k .eq.  knuc) then
          if(aTrack%p%subcode  .eq. regptcl) then
             if( aTrack%p%charge == 0 ) then
                ok = ke .gt. EminH
             else
                ok = ke > KEmin
             endif
          else
!              anti particle
             ok=aTrack%p%fm%p(4) > KEmin
          endif
       elseif( k .eq.  kpion .or.
     *         k .eq.  kkaon .or.
     *         k == kmuon ) then
!            can decay
          ok = aTrack%p%fm%p(4) > KEmin
       else
          if(k == kgnuc ) then
             ok = ke  > KEmin
          else
             ok = ke  >  KEmin
          endif
       endif
       if(ok) then
          icon=0
       else
          icon=1
          if( k == knuc .and. aTrack%p%charge == 0 ) then
             icon = -1
          endif
       endif
       end

      subroutine epAbsorb(aTrack, icon)
      use epModify
      use moddedx
      use modUI
      use modGUI
      implicit none
#include "ZmediaLoft.h"          
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
#include "Zcode.h"
#include "Zmass.h"

       type(epTrack)::  aTrack   ! input. this track's energy is examined
      integer icon              ! output. 0. the particle is still alive
                                !         1. death. 

      logical ok, kbtest, needdedx
      integer k
      real*8 ke, qf
      type(gui)::sinfo
      integer::modif
      real*8 cf, qg, cff
      integer,save::nstrange=0
!////////////////                                                              
!      logical show
!      common /showshow/ show
!///////////       
      qf = 0.
      ke = aTrack%p%fm%p(4)- aTrack%p%mass
      k =  aTrack%p%code
      icon = 1

      Move%dE= 0.
      Move%dEeff= 0.
      Move%dEioni = 0.

      if(k .eq. kelec) then
          if(aTrack%p%charge .eq. -1) then
             if(kbtest(Eabsorb, BitElectron)) then
                Move%dE = ke
                Move%dEeff = ke
                Move%dEioni= Move%dE  
                SumDe = SumDe + Move%dE
                call epLightPreUserde(1, aTrack)
             endif
          else
!              positron
             if(kbtest(Eabsorb, BitPositron)) then
                Move%dE = aTrack%p%fm%p(4) + masele
                SumDe = SumDe + Move%dE
                Move%dEeff = Move%dE
                Move%dEioni= Move%dE  
                call epLightPreUserde(1, aTrack)
             endif
          endif
       elseif(k .eq. kphoton) then
          if(kbtest(Eabsorb, BitPhoton)) then
             Move%dE = aTrack%p%fm%p(4)
             Move%dEeff = ke
             Move%dEioni= Move%dE    ! photo-electric electron's
             SumDe = SumDe + Move%dE
             call epLightPreUserde(1, aTrack)
          endif
       elseif(k .eq.  knuc) then
          if(aTrack%p%subcode  .eq. regptcl) then
             if( (aTrack%p%charge .eq. 1 .and. 
     *            kbtest(Eabsorb, BitProton)) .or.
     *            (aTrack%p%charge .eq. 0 .and.
     *            kbtest(Eabsorb, BitNeutron))) then
                Move%dE = ke
                Move%dEeff = ke
                Move%dEioni= Move%dE  
                SumDe = SumDe + Move%dE
                call epLightPreUserde(1, aTrack)
             endif
          else
!              anti particle
!             energy to be liberated further cannot be estimated.
!               (dependent on the ptcl type)
             if(kbtest(Eabsorb, BitAntiNuc)) then
                Move%dE = aTrack%p%fm%p(4) + aTrack%p%mass
                Move%dEeff = ke
                Move%dEioni= Move%dE  
                SumDe = SumDe + Move%dE
                call epLightPreUserde(1, aTrack)
             endif
          endif
       elseif( k .eq.  kpion .or.
     *         k .eq.  kkaon .or.
     *         k == kmuon ) then
!            can decay
          if(kbtest(Eabsorb,BitDecay)) then
             Move%dE = ke
             Move%dEeff = ke
             Move%dEioni= Move%dE  
             SumDe = SumDe + Move%dE
             call epLightPreUserde(1, aTrack)
          endif
       else
          if(k == kgnuc ) then
             if(   kbtest(Eabsorb, BitProton) ) then
                if( aTrack%p%charge > 1  ) then
                   needdedx = Media(MediaNo)%Birks /= ' '
                   modif = Det%cmp(Cn)%modifier
                   sinfo%modif = modif
                   if(.not. needdedx) then
                      if(modif > 0 .and.
     *                   allocated( modify) ) then
                         needdedx =
     *                        IBITS(modify( modif )%kind, bitQuench,
     *                          1)  > 0  
                      endif
                      if(needdedx ) then
!                        to get queching effect we need dedx
                         call epdedxhvy(Media(MediaNo), 
     *                   aTrack%p, dedx, dedxf)
                         if(HowQuench == 0 ) then
                            call epOrgCorrec(modif, 
     *                        Media(MediaNo), aTrack%p, dedx, cf)
                         elseif( HowQuench >= 1 ) then
                            call epOrgCorrec( modif,
     *                    Media(MediaNo), cTrack%p, dedxf, cff)
                            sinfo%dedx = dedx
                            sinfo%dedxf = dedxf
                            call epGUI(1, sinfo)
                            cf =  sinfo%cf
                            qf =  sinfo%qf
                         endif
                      else
                         cf =1.0
                      endif
                   else
                      cf = 1.0
                   endif
                   if( aTrack%p%charge /= 0  ) then  ! for safety
                      Move%dE = ke
                      Move%dEeff = ke *(cf*qf+ 1.-qf)
                      Move%dEioni= Move%dE  
                      SumDe = SumDe + Move%dE
                      call epLightPreUserde(1, aTrack)
                   endif
                endif
             endif
          else
             if(k /= kneue .and. k/= kneumu ) then
                if(kbtest(Eabsorb, BitOther)) then
                   Move%dE = ke
                   Move%dEeff = ke
                   Move%dEioni= Move%dE  
                   SumDe = SumDe + Move%dE
                   call epLightPreUserde(1, aTrack)
                endif
             endif
          endif
       endif


       if(Move%Abort .ne. 0) then
          if( Move%Abort == 3 ) then
             Move%Abort = 0
          else
             call epempty       ! empty the stack
             call epSkipUpdateNo
          endif
          icon = 1
       endif
       end
      subroutine epchckE(aTrack, icon)
      use epModify
      use moddedx
      use modGUI
      implicit none
#include "ZmediaLoft.h"          
#include "ZepTrackv.h"
#include "ZepTrackp.h"
#include "Zcnfig.h"
#include "Zcode.h"
#include "Zmass.h"
!            aTrack is examined if it's  energy is too low
!
       type(epTrack)::  aTrack   ! input. this track's energy is examined
      integer icon              ! output. 0. the particle is still alive
                                !         1. death. 

      logical ok, kbtest, needdedx
      integer k
      real*8 ke
      real*8 cf,qf,cff
      type(gui):: sinfo
      integer:: modif
!////////////////                                                              
!      logical show
!      common /showshow/ show
!///////////       

      ke = aTrack%p%fm%p(4)- aTrack%p%mass
      k =  aTrack%p%code
! >>>>>>>>>>>>>>>>>>light
      if(k .eq. klight) then
         if(Light <=  0 ) then
            write(0,*) ' Light =0 but light appeared'
            icon = 1   ! light is not treated, so death
            stop
         endif
         call epLightchkE( aTrack, icon) ! check wave length
         return  ! *************
      endif
!<<<<<<<<<<<<<<<<<<<<light
      if(k == kgnuc ) then
!         ke = ke/aTrack.p.subcode  ! dE is by total K.E
      endif
      
      if(k .eq. kelec) then
          if(aTrack%p%charge .eq. -1) then
             ok=aTrack%p%fm%p(4) .gt. EminElec
             if(.not. ok .and. ke .gt. 0.) then
                if(kbtest(Eabsorb, BitElectron)) then
                   Move%dE = ke
                   Move%dEeff = ke
                   Move%dEioni= Move%dE  
                   SumDe = SumDe + Move%dE
!                   if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>light
                      call epLightPreUserde(1, aTrack)
!                   endif                               <<<<<<<<<<
                endif
             endif
          else
!              positron
             ok=aTrack%p%fm%p(4) .gt. KEmin
             if(.not. ok) then
                if(kbtest(Eabsorb, BitPositron)) then
                   Move%dE = aTrack%p%fm%p(4) + masele

                   SumDe = SumDe + Move%dE
                   Move%dEeff = Move%dE
                   Move%dEioni= Move%dE  
!                   if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>>light
                      call epLightPreUserde(1, aTrack)
!                   endif                               <<<<<<<<<<< 

                endif
             else
                if(aTrack%p%fm%p(4) .lt. masele*1.001d0 ) then
!                   if energy is very low, forced anihilation 
                   Move%proc='anih'
                   Move%Trunc=.false.
                endif
             endif
          endif
       elseif(k .eq. kphoton) then
          ok= aTrack%p%fm%p(4) .gt. EminGamma
          if(.not. ok) then
             if(kbtest(Eabsorb, BitPhoton)) then
                Move%dE = aTrack%p%fm%p(4)
                Move%dEeff = ke
                Move%dEioni= Move%dE  
                SumDe = SumDe + Move%dE
!                if(Det.cmp(Cn).CountDE .ge. 1) then  >>>>>>>>>>>light
                   call epLightPreUserde(1, aTrack)
!                endif                                <<<<<<<<<<<
             endif
          endif
       elseif(k .eq.  knuc) then
          if(aTrack%p%subcode  .eq. regptcl) then
!             ok= ke .gt. KEmin
             if( aTrack%p%charge == 0 ) then
                ok = ke .gt. EminH
             else
                ok = ke > KEmin
             endif
             if(.not. ok) then
                if( (aTrack%p%charge .eq. 1 .and. 
     *             kbtest(Eabsorb, BitProton)) .or.
     *              (aTrack%p%charge .eq. 0 .and.
     *             kbtest(Eabsorb, BitNeutron))) then
                   Move%dE = ke
                   Move%dEeff = ke
                   Move%dEioni= Move%dE  
                   SumDe = SumDe + Move%dE
!                   if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>>>>>>>light
                      call epLightPreUserde(1, aTrack)
!                   endif                               <<<<<<<<<<<  
                endif
             endif
          else
!              anti particle
             ok=aTrack%p%fm%p(4) .gt. KEmin
             if(.not. ok) then
!                energy to be liberated further cannot be estimated.
!                (dependent on the ptcl type)
                if(kbtest(Eabsorb, BitAntiNuc)) then
                   Move%dE = aTrack%p%fm%p(4) + aTrack%p%mass
                   Move%dEeff = ke
                   Move%dEioni= Move%dE  
                   SumDe = SumDe + Move%dE
!                   if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>>>>>>>light
                      call epLightPreUserde(1, aTrack)
!                   endif                               <<<<<<<<<<<<<
                endif
             endif
          endif
       elseif( k .eq.  kpion .or.
     *         k .eq.  kkaon .or.
     *         k == kmuon ) then
!            can decay
!          ok = aTrack.p.fm.p(4) .gt. KEmin
          ok = aTrack%p%fm%p(4) > EminH
          if( ok .and. ke <= EminH ) then
             if( k == kmuon .and.  aTrack%p%charge == -1) then
                ! can be abosrbed by capture so follow until death
             elseif( .not. Move%trunc ) then
                ! force to decay
                Move%proc='decay'
                call cresetIntInf   ! this rest is needed
             endif
          endif
          if(.not. ok) then
             if(kbtest(Eabsorb,BitDecay)) then
                Move%dE = ke
                Move%dEeff = ke
                Move%dEioni= Move%dE  
                SumDe = SumDe + Move%dE
!                if(Det.cmp(Cn).CountDE .ge. 1) then  >>>>>>>>>>>>light 
                   call epLightPreUserde(1, aTrack)
!                endif                                <<<<<<<<<<<<<<<
             endif
          endif
       else
          if(k == kgnuc ) then
             ok = ke  > KEmin
             if(.not. ok) then
                if(   kbtest(Eabsorb, BitProton) ) then
                   if( aTrack%p%charge > 1  ) then
                      needdedx = Media(MediaNo)%Birks /= ' '
                      modif = Det%cmp(Cn)%modifier
                      sinfo%modif = modif
                      if(.not. needdedx) then
                         if(modif > 0 .and.
     *                      allocated( modify) ) then
                            needdedx =
     *                      IBITS(modify(modif )%kind,
     *                      bitQuench, 1) > 0  
                         endif
                      endif
                      if(needdedx ) then
!                        to get queching effect we need dedx
                         call epdedxhvy(Media(MediaNo), 
     *                   aTrack%p, dedx, dedxf)
                         if( HowQuench == 0 ) then
                            call epOrgCorrec(modif, 
     *                       Media(MediaNo), aTrack%p, dedx, cf)
                         elseif( HowQuench >= 1 ) then
!                                next call is not needed now
                            call epOrgCorrec( modif,
     *                    Media(MediaNo), cTrack%p, dedxf, cff)
                            sinfo%dedx = dedx
                            sinfo%dedxf = dedxf
                            call epGUI(1, sinfo)
                            cf = sinfo%cf
                            qf = sinfo%qf
                         endif
                      else
                         cf =1.0
                      endif
                   else
                      cf = 1.0
                   endif
                   if( aTrack%p%charge /= 0  ) then  ! for safety
                      Move%dE = ke
                      Move%dEeff = ke *(cf*qf+1.0-qf)
                      Move%dEioni= Move%dE  
                      SumDe = SumDe + Move%dE
                      call epLightPreUserde(1, aTrack)
                   endif
                endif
             endif
          else
             ok = ke  >  KEmin
             if(.not. ok) then 
                if(k /= kneue .and. k/= kneumu ) then
                   if(kbtest(Eabsorb, BitOther)) then
                      Move%dE = ke
                      Move%dEeff = ke
                      Move%dEioni= Move%dE  
                      SumDe = SumDe + Move%dE
!                if(Det.cmp(Cn).CountDE .ge. 1) then >>>>>>>>>>>>light
                      call epLightPreUserde(1, aTrack)
!                endif                               <<<<<<<<<<<<
                   endif
                endif
             endif
          endif
       endif
       if(ok) then
          icon=0
       else
          icon=1
       endif
       if(Move%Abort .ne. 0) then
          if( Move%Abort /= 3 ) then
             call epempty       ! empty the stack
             call epSkipUpdateNo
          else
             Move%Abort = 0
          endif
          icon = 1
       endif
       end
!      *****************
       subroutine epaddTime
!           >>>>>>>>>>>>>>>light
       use modepLightPty
!           <<<<<<<<<<<<<<<
       implicit none
!        new position is assumed to be fixed.
!        (however, before scattering, and energy check).
!        update time and take trace
#include  "Zglobalc.h"
#include  "ZepTrackp.h"
#include  "Zmedia.h"           
#include  "ZepTrackv.h"

       real*8 ctau, u
       real*8 beta1, betaav
!
       if(Move%Track%p%mass .eq. 0.) then
!   >>>>>>>>>>>>>>>>>>>light
          if(Move%Track%p%code < 0 .and. cLcompNo > 0 ) then
             betaav = 1./Lcomp( cLcompNo )%refracN
          else
!  <<<<<<<<<<<<<<<<<<<<<<<
             betaav=1.d0
          endif
       else
          call cgetBeta( Move%Track%p,  betaav)
          if(betaav .lt. 0.98) then
             call cgetBeta(cTrack%p, beta1)
             betaav = (betaav+beta1)/2.
          endif
       endif
       if(betaav .gt. 0.) then
          Move%Track%t = cTrack%t + Move%dl/betaav
       else
!               stopped one. add decay time if possible
          call cgetctau(cTrack%p, ctau)
          if(ctau .eq. Infty) then
!               If the particle is still to be treated
!               even after it stops, the ptcl 
!               should be antiprpton or anti-somthing to
!               anihilate in the next step. We assume
!               the anihilation takes place instantly
!               after stopping. so don't add any time
!               
!                Move.Track.t= cTrack.t+1.e8  ! version <= 8.62
                Move%Track%t = cTrack%t
          else
             call rndc(u)
             Move%Track%t = cTrack%t - log(u)*ctau
          endif
       endif
       end

      subroutine epl2wTrack(aTrack, bTrack)
!        convert local info (pos, w, p) into world ones
      implicit none
#include  "ZepTrack.h"

       type(epTrack)::  aTrack ! input
       type(epTrack)::  bTrack ! output.
      bTrack = aTrack
      call epl2w(aTrack%cn, aTrack%pos, bTrack%pos)      
      call epl2wdm(aTrack%cn, aTrack%w, bTrack%w, bTrack%p)
      end    subroutine epl2wTrack
!         returns dE/dx (GeV/cm2)  ; this should not be moved
!       to epquery (due to use moddedx)
      subroutine epqElossRate(dedxout)
      use moddedx
      implicit none
      real(8),intent(out):: dedxout
      dedxout = dedx
      end

      subroutine epe2p(aTrack)
!         E & direc cos  -->  px,py,pz
      implicit none
#include "ZepTrack.h"
      
       type(epTrack)::  aTrack

      real*8 p

      p = aTrack%p%fm%p(4)**2 - aTrack%p%mass**2
      if(p .lt. 0.) then
         p = 0.
      else
         p = sqrt(p)
      endif

      aTrack%p%fm%p(1) = p * aTrack%w%x
      aTrack%p%fm%p(2) = p * aTrack%w%y
      aTrack%p%fm%p(3) = p * aTrack%w%z
      end
      subroutine epResetCountIO(cmpNo, countio) 
      implicit none
#include  "Zmedia.h"          
#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in):: cmpNo !  comp. #                                                                      
      integer,intent(in):: countio  ! countio to be set                                                           
      if(cmpNo >= 1 .and. cmpNo <= Det%nct) then
         Det%cmp(cmpNo)%CountIO = countio
      else
         write(0,*) 'Warning from epResetCountIO:'
         write(0,*) ' specified comp.# ', cmpNo, ' non exsistent'
      endif

      end
#if !defined (INTINFO)
!          this is dummy routine to avoid link problem
!      subroutine epUI(info, loc1, loc2)
!      implicit none
!      integer,intent(in):: info
!      integer,intent(in):: loc1,loc2
!      end
#endif

