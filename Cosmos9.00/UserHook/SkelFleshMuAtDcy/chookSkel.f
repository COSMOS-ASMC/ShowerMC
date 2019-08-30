!   copy this as chookSkel.f; make -f chookSkel.mk
!    This gives
!    1)  Emu of observed muons
!    2)  eta of parent( pi or K)
!    3)  Xcm of parent (pi or K)
!    4)  index to show the parent is  pi or K.
!    5)  the number of collisions the ancestors suffered.
! To do so, at each interaction, we calculate eta, Xcm of
!    produced pi or K.  Count number of collisions
!    and embed these and index in .user; For observed muons,  we extract
!    these quantities.    
!
#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
!
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      implicit none
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Zprivate.h"

      character*10  append
      character*100 msg
      logical  ex, opn
      integer klena

!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)


      call cqUHooki(1, Mdev)      ! get skeleton memo dev #
      call cqUHooki(2, Wdev)      ! get working disk dev #


      call cqUHooki(3, NgMin)     ! get Nh min
      call cqUHooki(4, NhMin)     ! get Ng min
      call cqUHooki(5, Where)     ! where to check 

      call cqUHookr(1, SumegMin)      ! sum E min
      call cqUHookr(2, SumehMin)

      

      call cqUHookc(1, msg)       ! get file name for sekelton memo
      call cgetfname(msg, Mskel)  ! add host name etc if needed
      call cqUHookc(2, msg)       ! get file name for working
      call cgetfname(msg, Wskel)  ! add host name etc if neededn
      call cqUHookc(3, append)    ! append data, if Mskel already exists

      write(msg, *) 'Skeleton is  judged at obs.pos=', Where
      call cerrorMsg(msg, 1)
      write(msg, *) ' Ngmin=',NgMin,
     *     ' SumEGmin=',SumegMin/1000.,'TeV'
      call cerrorMsg(msg, 1)
      write(msg, *) ' Nhmin=',NhMin,
     *    ' SumEHmin=',SumehMin/1000.,'TeV'
      call cerrorMsg(msg, 1)

!     
      inquire(file=Mskel(1:klena(Mskel)), opened=opn, exist=ex)
      if(opn) then
         call cerrorMsg(Mskel, 1)
         call cerrorMsg(' already opened: starange', 0)
      elseif(ex .and.
     *      append(1:klena(append)) .eq. 'append' ) then
         open(Mdev, file=Mskel, form='unformatted',status='old')
         call cerrorMsg('skeleton node info. will be appended', 1)
!           skip to the end of file
         do while( .true. )
            read(Mdev, end=100)
         enddo
      else
         if(ex .and.
     *      append(1:klena(append)) .ne. 'append' ) then
            call cerrorMsg(
     *           'Old skeleton node info. file exists', 1)
            call cerrorMsg(
     *           'but node info. will NOT be appended', 1)
         endif
         open(Mdev, file=Mskel(1:klena(Mskel)), form='unformatted',
     *       status='unknown')
      endif

 100  continue

      open(Wdev, file=Wskel(1:klena(Wskel)), form='unformatted',
     *      status='unknown' )

      Accepted = 0   ! counter;  accepted as skeleton 
!
      end


!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      implicit none

#include "Ztrack.h"
#include "Ztrackv.h"
#include "Ztrackp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include "Zcode.h"
#include "Zprivate.h"      



      type(coord )::angle


      integer  EventNo
      integer i, j
      integer seed(2)
      real*8 svEasWait, svEthin(2), kepn

      bgevent = .true.
      Np = 0
      call cqIncident( incident, angle)
      kepn = incident.p.fm.p(4)
      if( incident.p.code .eq. kgnuc ) then
         kepn = kepn/ incident.p.subcode
      endif
      Ethresh = kepn * WaitRatio

      svEasWait = EasWait       ! for safety save
      svEthin(1) = Ethin(1)           ! //
      svEthin(2) = Ethin(2)           ! //
      call csetEmin(Generate2, KEminObs2, Cutneg, Cuteg)
      EasWait = svEasWait       ! restore
      Ethin(1) = svEthin(1)
      Ethin(2) = svEthin(2)

      rewind  Wdev

!     ===================================================

      EventNo = EventNo + 1     
      do i = 1, NoOfSites
         write(*, 1000)
     *    sngl(ObsSites(i).pos.depth),
     *    EventNo, 
     *    incident.p.code,
     *    incident.p.subcode,
     *    incident.p.charge, 
     *    incident.p.fm.e,
     *    -angle.r(1),
     *    -angle.r(2),
     *    -angle.r(3)
 1000       format(f10.3,i9,3i4,e15.5,3(1x,f12.8))
      enddo


!     ===================================================
      end
  

!     ************************************ hook for observation
!     *  One particel information is brought here by the system.
!     *  All information of the particle is in aTrack
!     *
      subroutine chookObs(aTrack, id)
!
!     Note that every real variable is in double  precision so
!     that you may output it in sigle precision to save the memory.
!     In some cases it is essential to put it in sigle (say,
!     for gnuplot).
! 
      implicit none
#include "Zcode.h"
#include "Ztrack.h"
#include "Zprivate.h"
      integer id  ! input.  2 ==> reached at an observation level
!                           1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
!     ===================================================
      common /testcos/sumg, ng(20), nth, EventNo
      real*8 sumg
      integer ng, nth, EventNo

!     ===================================================
      real*8 user
      real*4 userE(2)
      integer*2 userI(4)
      equivalence (user, userE(1)), (user, userI(1))

      integer codex
      logical kbitest
      save
                        
      if( id .eq. 2 .and. aTrack.p.code .ne. kneumu .and.
     *   aTrack.p.code .ne. kneue) then
         Np = Np + 1
         if( Np .gt. NpMax) then
            call cerrorMsg(
     *   '# of particles >NpMax in observation', 0)
         endif
         o(Np).where = aTrack.where
         o(Np).code = aTrack.p.code
         o(Np).subcode = aTrack.p.subcode
         o(Np).charge = aTrack.p.charge
         o(Np).atime = aTrack.t
         o(Np).erg = aTrack.p.fm.p(4)
         o(Np).mass = aTrack.p.mass
         o(Np).x = aTrack.pos.xyz.r(1)
         o(Np).y = aTrack.pos.xyz.r(2)
         o(Np).wx = aTrack.vec.w.r(1)
         o(Np).wy = aTrack.vec.w.r(2)
         o(Np).wz = aTrack.vec.w.r(3)
         o(Np).zenith = aTrack.vec.coszenith
         o(Np).user =  aTrack.user
!     ===================================================
!       if( o(Np).code .le. 6 .and. o(Np).code .ne. 3 ) then
         if(o(Np).code .eq. 3) then
            user=aTrack.user
            call getpcode(userE(2), codex)
            write(*, 959) 
!     *  o(Np).where,   
!     *  o(Np).code,   
!     *  o(Np).charge,  
     *  ( o(Np).erg - o(Np).mass ),   
!     *  ( o(Np).x ), ( o(Np).y ) ,  
!     *  ( o(Np).wx ),  
!     *  ( o(Np).wy ),  
!     *  ( o(Np).wz ),  
!     *  ( o(Np).zenith ),
     *  userI(1)/1000., userI(2), userE(2), codex

! 959    format(3i3, g14.3,2f16.6,4(1x,f12.8), i3, 1pE11.3)
 959    format(2g14.3,i7, E11.3, i3)
      endif

!     ===================================================
      endif
!     ===================================================

!     ===================================================
      end

!    *********************************** hook for end of 1 vent
!    * At this moment, 1 event generation has been ended.
!    *
      subroutine chookEnEvent

      implicit none
#include "Ztrack.h"
#include "Zcode.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include "Zmanagerp.h"
#include "Zprivate.h"

      integer i
      integer seed(2)
      real  sumeg, sumeh
      logical memorize

      integer ng, nh
      
      ng = 0
      nh = 0

      sumeg = 0
      sumeh =0
!         count sum Eg etc
      do i = 1, Np
         if(o(i).where .eq. Where) then
            if(o(i).code .le. kelec) then
               ng = ng + 1
               sumeg = sumeg + o(i).erg
            elseif( ( o(i).code .ge. kpion .and.
     *                o(i).code .le. knuc )
     *             .or.
     *               o(i).code .eq. kgnuc) then
               nh = nh + 1
               sumeh = sumeh + o(i).erg
            endif
         endif
      enddo

!     ===================================================
      memorize =(ng .ge. NgMin .and. sumeg .ge. SumegMin) .or.
     *  ( nh .ge. NhMin .and. sumeh .ge. SumehMin)

!/////
      write(0,*) ' memo=', memorize,
     *   ' ng=',ng, ' seg=',sumeg, ' nh=',nh,' seh=',sumeh
!//////


!     ===================================================
      if(memorize) then
         Accepted = Accepted + 1
         call cwriteSeed        !  SeedFile
!          flag for end of 1 event on working disk
         write(Wdev)  -1, p
         call cmemorize(Wdev, Mdev)     !  reocord this event
      endif
      end


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun
      implicit none
#include "Ztrack.h"
#include "Zprivate.h"
      character*100 msg
!     =========================================


!     =========================================
      call cerrorMsg('++++++++++++', 1)
      write(msg, '(i8, a)') Accepted,
     *       ' events are memorized as skeleton'
      call cerrorMsg(msg, 1)
      call cerrorMsg('their  seeds are also memorized', 1)
      call cerrorMsg('-----------',1)
      call  cprintStatus   ! if don't like,  comment out
      end
!     ********************************* hook for trace
!     *  This is called only when trace > 100
!     *  User should manage the trace information here.
!     *  If you use this, you may need some output for trace
!     *  at the beginning of 1 event generatio and at the end of  1 event
!     *  generation so that you can identfy each event.
!     *
!     *
      subroutine chookTrace
            implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsv.h"

       real*4 h1,  h2
!
!    Every time a particle is moved in the atmosphere, this routine is called,
!    if trace > 100
!         For a one track segment,
!     TrackBefMove  has  track information at the beginning of the segment.
!     MoveTrack    has   track information at the end of the segment.
!   
!     You can know the  information a track contains in the 
!     chookObs routine. (Note however, no conversion of coordinate
!     has been done.  The values are in the Earth xyz system.)
!     Besides quantities explained there, you can use, for a  given 'track'
!
!     atrack.pos.xyz.x, atrack.pos.xyz.y, atrack.pos.xyz.z    (x,y.z)
!     atrack.pos.radiallen   (distance from the center of the earth)
!     atrack.pos.depth       (vertical depth)
!     atrack.pos.height      (vertical heigth from sea level)  
!

      h1 = TrackBefMove.pos.height- ObsSites(NoOfSites).pos.height
      h2 = MovedTrack.pos.height - ObsSites(NoOfSites).pos.height

      end
!     ********************* this is the hook called when
!       an electron made an interaction.
!
      subroutine chookEInt(never)
      implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zcode.h"
#include  "Zprivate.h"
!  #include  "Ztrackp.h"
      integer never             ! input & output
      integer nlow
! 
!    If Job = 'newskel', input  "never" may be < 0,  and MovedTrack
!    may not be an electron.
!      never = -1 :  MovedTrack K.E becomes < KEminObs due to
!                    energy loss during traveling.
!      never = -2 :  The same as above, but the particle crosses an
!                    observation depth, so the calling point to this
!                    routine is different from the never = -1 case.
!      never = -3 :  K.E >= KEminiObs.  The ptcl is observed at the
!                    current deepest Obs. level. but at the flesh time
!                    the deepest level will be more deep so that
!                    this must be memorized.
!
!         For above cases, the product is only MovedTrack and 
!         no particle is in PWork.
!  Otherwise,
!   MovedTrack is the electron that made interaction
!   Pwork contains produced particles (normally gamma, but mayb  e).
!   Nproduced has the number of particles in Pwork (normally 1)
!   IntInfArray(ProcessNo) contains the type of interaction
!        IntInfArray(ProcessNo).process will have one of
!       'brems', 'mscat', 'bscat'  'anihi' or 'mbrem'
!     
      if(never .lt. 0) then
         Nproduced = 1
!           IBM user must modify next
         Pwork(1) = MovedTrack.p
      endif

!         high energy parent at node might be used
!        for hybrid AS generation if it is in some
!        energy region.
      if( MovedTrack.p.code .eq. kelec ) then
         if( MovedTrack.asflag .eq. 0 ) then
            if( MovedTrack.p.fm.p(4) .lt. Ethresh ) then
               MovedTrack.asflag = -1
            endif
         endif
      endif

      call cmemoNode(Wdev, never, 2)
      if(MovedTrack.asflag .eq. -1) then
!          decendent should not be used to  generae A.S
         MovedTrack.asflag  = -2
      endif

 10   continue
      never = 0
      end

!     ********************* this is the hook called when
!       a gamma ray made an interaction.
!
      subroutine chookGInt(never)
      implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zprivate.h"
!  #include  "Ztrackp.h"
      
      integer never   ! input & output
      
!         don't make never = 1, if you want to get
!         information after a gamma ray made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the gamma that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!         IntInfArray(ProcessNo).process will have one of
!        'pair', 'comp', 'photoe' 'photop' 'mpair'
      call cmemoNode(Wdev, 1, 1)
      never = 0
      end

!     ********************* this is the hook called when
!       non e-g particle made an interaction.
!
      subroutine chookNEPInt(never)
      implicit none
#include  "Zmaxdef.h"
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zstackv.h"
#include  "Zprivate.h"
      
      integer never   ! input & output

      type(coord):: angle
      type(ptcl)::  aptcl
      type(ptcl):: tgt, pj
      real*8 user
      real*4 userE(2)
      integer*2 userI(4)
      equivalence (user, userE(1)), (user, userI(1))
      real*8 y, eta
      integer i, icon, codex
      integer stackpos
      save 
      
      if(bgevent) then
         user=MovedTrack.user
         userI(2) = 0
         userE(2) = 0.
         MovedTrack.user = user
         bgevent = .false.
      endif

      proc = IntInfArray(ProcessNo).process

      incip = MovedTrack.p
      if(proc .eq. "coll") then
         call cmkptc(6, 0, 1, tgt)
         tgt.fm.p(1) = 0.
         tgt.fm.p(2) = 0.
         tgt.fm.p(3) = 0.
         tgt.fm.p(4) = tgt.mass
         pj = incip
         pjcm = pj
         call cirot3vec(1, incip, incip, pj)
         call cgeqm(pj, tgt, cmsptcl, icon)
         call cbst1(1, cmsptcl, pj, pjcm)
      endif

      call cgetCurrentStackPos(stackpos)
      do i = stackpos-Nstacked+1, stackpos
         user = Stack(i).user
         if(proc .eq. "coll") then
            aptcl=Stack(i).p
            call cirot3vec(1, incip, aptcl, aptcl)
            call cbst1(1, cmsptcl, aptcl, aptcl)
!                y, eta
            call cyeta(aptcl, y, eta)
            userI(2)=userI(2)+1
            userI(1)=eta*1000.
!                xcm;  if anti-proton or anti-neutron
!                it may be stopped one
            if(pjcm.fm.p(3) .eq. 0.) then
               userE(2)= aptcl.fm.p(3)/(cmsptcl.mass/2.)
            else
               userE(2)= aptcl.fm.p(3)/pjcm.fm.p(3)
            endif
            Stack(i).user=user 
         elseif( proc .eq. "decay"  .and.
     *      ( incip.code .eq. 4 .or. incip.code .eq. 5))  then
            call setpcode(incip, userE(2))
            Stack(i).user=user 
         endif
      enddo
!         don't make never = 1, if you want to get
!         information after a non-e-g particle  made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the particle that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
      if( proc .eq. "coll"  ) then
         call cmemoNode(Wdev, 1, 4)
      else
         call cmemoNode(Wdev, 1, 3)
      endif

      end
      subroutine cmemoNode(dev, flag, neporeg)
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zprivate.h"
!
      integer dev  ! input. fortran logical dev. number for data write
      integer flag ! input. If this is -3, high energy particles must
                   !       also be memorized.
      integer neporeg ! input 1; from g 2 from e; 3 from NEP
!
!
!       memorize nodal information at interaction.
!     
!        \
!         \  projectile    = MovedTrack
!         /|\  produced particles.: Pwork(i), i=1~Nproduced
!        / | \
!
!   From projectile, the track information is extracted
!   From produced particles, only those of K.E< KEminObs is
!   extracted and information needed for further tracking is 
!   obtaned and memorized. The position information is common
!   to all the children.
!
!   Track   information;   pos =
!                                xyz = r(1~3), sys
!                                radiallen, depth, height, colheight
!                            t
!                          vec  =
!                                 w = r(1~3), sys
!                                 coszenith
!                          wgt
!                         where
!                        asflag 
!
!      Among these, we don't memorize 
!         sys which is always 'xyz'
!       radiallen: can be computed from height
!         vec; children knows  their direction
!         wgt: normally not needed, it should be 1 for skeleton making
!              So thinsamping must not be used when making skeleton.
!      asflag: should be always F, (assume for skeleton making, hybrid
!              air shower is not requested)       
!
!   write  track info. of projectile
!
      integer i, nlow
      real*8 ke
      real*8  y, eta
      real*8 user
      real*4 userE(2)
      integer*2 userI(4)
      equivalence (user, userE(1)), (user, userI(1))
      type(ptcl):: aptcl
      integer codex
      save

!
!
      nlow = 0 
      do i = 1, Nproduced
         ke = Pwork(i).fm.p(4) - Pwork(i).mass 
         if( ( Pwork(i).code .le. kelec .and. 
     *       ke .ge.  Cuteg ) .or. 
     *       ( Pwork(i).code .gt. kelec .and.          
     *       ke .ge.  Cutneg ) ) then
!              count low energy ptcls
            if(flag .ne. -3) then
               if( ke .lt.  KEminObs) then
                  nlow = nlow + 1
               endif
            else
!               all ptcl must be memorized
               nlow = nlow + 1
            endif
         endif
      enddo
!         next is harmful for B-approx air shower
!         because, parent is not memorized and cannot
!         generate AS at flesh time at proper positon
!         asflag == -1 must apear in flesh
!      if(nlow .eq. 0 )  return !  *************
!       if there is no return here, the size of skeleton data
!       becomes 3.2 times larger

      if(nlow .eq. 0  .and. MovedTrack.asflag .ne. -1) return


      p.posx = MovedTrack.pos.xyz.r(1)
      p.posy = MovedTrack.pos.xyz.r(2)
      p.posz = MovedTrack.pos.xyz.r(3)
      p.depth = MovedTrack.pos.depth
      p.height = MovedTrack.pos.height
      if( MovedTrack.pos.colheight .gt. 1.e36) then
         p.colHeight = 1.e36    ! no col. yet.
      else
         p.colHeight =  MovedTrack.pos.colheight
      endif
      p.atime = MovedTrack.t
      p.where = MovedTrack.where
      p.coszenith = MovedTrack.vec.coszenith
      p.code = MovedTrack.p.code
      p.erg  = MovedTrack.p.fm.p(4)
      p.asflag = MovedTrack.asflag
      p.user =  MovedTrack.user

!c
!c     *           p.posx, p.posy, p.posz, p.depth, p.height, 
!c     *           p.colHeight, p.atime, p.where
!c
!       write particle info
!           p(1~4)
!           mass
!           code
!           subcode
!           charge
! 


      write(dev) nlow, p

      do i = 1, Nproduced
         ke = Pwork(i).fm.p(4) - Pwork(i).mass 
         if( ( Pwork(i).code .le. kelec .and. 
     *       ke .ge.  Cuteg ) .or. 
     *       ( Pwork(i).code .gt. kelec .and.          
     *       ke .ge.  Cutneg ) ) then
            if(flag .eq. -3 .or.  ke .lt. KEminObs) then
               c.code = Pwork(i).code
               c.subcode =  Pwork(i).subcode 
               c.charge =  Pwork(i).charge
               c.fm(1) = Pwork(i).fm.p(1)
               c.fm(2) = Pwork(i).fm.p(2)
               c.fm(3) = Pwork(i).fm.p(3)
               c.fm(4) = Pwork(i).fm.p(4)
               c.mass = Pwork(i).mass
               if(neporeg .eq. 4 ) then
                  aptcl = Pwork(i)
                  call cirot3vec(1, incip,  aptcl,  aptcl)
                  call cbst1(1, cmsptcl, aptcl, aptcl)
                  call cyeta(aptcl, y, eta)
                  user=MovedTrack.user
                  userI(1)=eta*1000.
                  userI(2) = userI(2) + 1
!                    xcm
                  if( pjcm.fm.p(3) .eq. 0.) then
                     userE(2)= aptcl.fm.p(3)/(cmsptcl.mass/2.)
                  else
                     userE(2)= aptcl.fm.p(3)/pjcm.fm.p(3)
                  endif
                  c.user = user
               elseif( proc .eq. "decay"   .and.
     *             incip.code .eq. 4 .or. incip.code .eq. 5) then
!                     set parent ptcl code by using last bit
                  user=MovedTrack.user
                  call setpcode( incip, userE(2) )
                  c.user = user
!////////////////////
!                  if(userE(2) .gt. -1. .and. 
!     *              Pwork(i).code .le. 5) then
!                     write(*,'(a,3i3, 4g11.3)')
!     *                 'l: ', Pwork(i).code, Pwork(i).subcode,
!     *                    Pwork(i).charge,
!     *                    eta, y, userE(2), Pwork(i).fm.p(4)
!                  endif
!/////////////////
               else
                  c.user = p.user
               endif
               write(dev) c
            endif
         endif
      enddo
      end

!        
      subroutine cmemorize(from, to)
      implicit none
#include  "Ztrack.h"
#include  "Ztrackv.h"

      integer from      !  working disk file 
      integer to        !  permanent disk file where skeleton is sotered




      integer num, cumnum, irevent(2)
!      type(track):: incident
!      type(coord):: angle


      rewind from
!          need not memorize, can be generated at flesh
!      call cqIncident(incident, angle)      

      call cqEventNo(num, cumnum)

      call cqIniRn(irevent)   ! seed of the event

      write(to) cumnum, num, irevent,  Zfirst
      call cputHES(to)   ! put high energy showers.
      write(0,*) ' first Z=',Zfirst.pos.depth*0.1, ' g/cm2',
     *  Zfirst.pos.height, ' m' 
      
      call cputNodInfo(from, to)  ! put nordal info.

      end

      subroutine cputHES(to)
      implicit none
#include "Ztrack.h"
#include "Zprivate.h"
      integer to
!
!        write high energy sekeleton data into 'to'      
!
      integer i

      write(to) Np

      do i = 1, Np
         write(to) o(i)
      enddo

      end

      subroutine cputNodInfo(from, to)
      implicit none
#include "Ztrack.h"
#include "Zprivate.h"

      integer from, to
      
      integer i, nlow

      nlow = 1
      do while ( nlow .ge. 0 )
         read(from) nlow,  p
         write(to)  nlow,  p
         do i = 1, nlow
            read(from) c
            write(to) c
         enddo
      enddo

      end


