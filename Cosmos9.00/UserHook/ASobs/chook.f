#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      implicit none
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"

!          ///////////
      integer nl, nth
      parameter (nl = 20, nth=12)
      common /testcos/eth(nth),
     * ng(nth, nl), ne(nth, nl), nmu(nth, nl), 
     * nh(nth, nl), ntha
      integer ng, ne, nmu, nh, ntha
      real*8 eth
!     //////////////

!         If you feel writing the parameters on stderr is
!         a bother, comment out the next or
!         use other device than ErrorOut.
!         Also you may comment out all output routines below.
#ifdef sun4
      external csigHandler
      integer  ieeer, ieee_handler

      ieeer = ieee_handler('set', 'invalid', csigHandler)
#endif

!
!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)

      eth(1) = 0.3d-3
      eth(2) = 0.5d-3
      eth(3)=  1.d-3
      eth(4) = 2.d-3
      eth(5) = 5.d-3
      eth(6) = 10.d-3
      eth(7) = 20.d-3
      eth(8)=  50.d-3
      eth(9) = 100.d-3
      eth(10) = 200.d-3
      eth(11) = 500.d-3
      eth(12) = 1.d0

!      ////////////
      ntha = 1  !  for each ptcle out put use only 1 threshold
!     /////////////
      end
#ifdef sun4
      integer function csigHandler(sig, code, context)
      implicit none
#include "Zmanagerp.h"
      integer sig, code, context(5)
      write(ErrorOut, *)  ' f.p exception content=' , context(4)
!      call abort()
      end
#endif

!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      implicit none
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"

!          ///////////
      integer nl, nth
      parameter (nl = 20, nth=12)
      common /testcos/eth(nth),
     * ng(nth, nl), ne(nth, nl), nmu(nth, nl), 
     * nh(nth, nl), ntha
      integer ng, ne, nmu, nh, ntha
      real*8 eth
!     //////////////
      type(track):: inci
      type(coord):: angle

      integer  nev
!     //////////////
      integer i, j
      integer seed(2)
      do i = 1, nl
         do j = 1, ntha
            ng(j, i) = 0
            ne(j, i) = 0
            nh(j, i) = 0
            nmu(j, i) = 0
         enddo
      enddo

!      write(*, *)
      call cqIncident(inci, angle)
      write(*,'(i7,i4,g13.4,3f10.7)') EventNo, inci.p.code, inci.p.fm.e,
     *  -angle.r(1),  -angle.r(2), -angle.r(3)      
      call cqIniRn(seed)
!      write(*,*) ' seed=', seed
      end
      subroutine ccount(nc, aTrack)
      implicit none
#include "Zcode.h"
#include "Ztrack.h"

      type(track):: aTrack
      integer i
!          ///////////
      integer nl, nth
      parameter (nl = 20, nth=12)
      common /testcos/eth(nth),
     * ng(nth, nl), ne(nth, nl), nmu(nth, nl),
     * nh(nth, nl), ntha
      integer ng, ne, nmu, nh, ntha
      real*8 eth
!
      integer nc(nth, nl)

      do i = 1, ntha
        if( aTrack.p.fm.e- aTrack.p.mass .lt. eth(i)) goto 10
        nc(i, aTrack.where) =  nc(i, aTrack.where) + 1
      enddo
 10   continue 
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
#include "Ztrackv.h"
#include  "Zheavyp.h"
      integer id  ! input.  1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack
      type(track):: inci
      type(coord):: angle, tetafai
!      integer i
!          ///////////

!          ///////////
      integer nl, nth
      parameter (nl = 20, nth=12)
      common /testcos/eth(nth),
     * ng(nth, nl), ne(nth, nl), nmu(nth, nl), 
     * nh(nth, nl), ntha
      integer ng, ne, nmu, nh, ntha
      real*8 eth
!     //////////////
      real*8 depth
      integer iij



!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
!      if(id .eq. 2 .and. aTrack.p.code .eq. kmuon ) then
       if( id .eq. 2) then
           call cqIncident(inci, angle)
           iij = aTrack.p.code  
!          call cgpid(iij, ptclid)
           if(iij .eq. kelec ) then
              call ccount(ne, aTrack)
!              ne(aTrack.where) = ne(aTrack.where) + 1
           elseif(iij  .eq. kphoton ) then
              call ccount(ng, aTrack)
!              ng(aTrack.where) = ng(aTrack.where) + 1
           elseif(iij .eq. kmuon ) then
              call ccount(nmu, aTrack)
!              nmu(aTrack.where) = nmu(aTrack.where) + 1
           elseif( iij .eq. kpion  .or. iij .eq.  kkaon .or. 
     *        iij .eq. knuc) then
              if(aTrack.p.charge .ne. 0 ) then
                 call ccount(nh, aTrack)
!                 nh(aTrack.where) = nh(aTrack.where) + 1
              endif
           endif
        endif

!            output typical quantities.
!         id .eq. 2 below  if want otuput
!        if(id .eq. 2 .and. aTrack.where .eq.  8 ) then
!        if(id .eq. 2 .and. aTrack.where .eq.  1 ) then
        if(id .eq. 0 ) then  ! never happens
!         if(aTrack.p.code .ne. kneue .and. aTrack.p.code .ne. 
!     *    kneumu) then
           call cwhere2dep(aTrack.where, depth)
!           write(*,*) int(depth/10. +0.001),
           write(*,*)
     *     aTrack.p.code,
!     *  aTrack.p.subcode,
!     *     aTrack.p.charge, 
!     *     sngl(aTrack.p.fm.e-aTrack.p.mass),
!     *      sngl(aTrack.vec.coszenith),
!     *      sngl(inci.p.fm.p(4)-inci.p.mass)/
!     *      Code2massN(int(inci.p.code)),
!     *      inci.p.code,
!     *     aTrack.pos.xyz.x, aTrack.pos.xyz.y,  aTrack.pos.xyz.z, 
     *    sngl(aTrack.vec.w.x), sngl(aTrack.vec.w.y), 
     *    sngl(aTrack.vec.w.z), 
     *     aTrack.vec.w.x**2+aTrack.vec.w.y**2+aTrack.vec.w.z**2
!     
!       write(*,
!     *  '(i2,1x,i2,1x,f12.2, g13.4, f12.2,1x, f12.2,1x,f7.4,i3)')
!     * nev, 
!     *  aTrack.where,   !  observation level. integer*2.  1 is highest.
!     *  aTrack.p.code,  ! " ", ptclid,    !  ptcl code.  integer*2.
!     *  aTrack.p.charge,  !  charge,  integer*2 
!     *  sngl(aTrack.t), !  relateive arrival time in nsec (NOT sec).
!                        !  if TimeStructure is F, nonsense.
!     *  sngl(aTrack.p.fm.e), !  - aTrack.p.mass), ! kinetic energy in GeV
!     *  sngl(aTrack.pos.xyz.x), sngl(aTrack.pos.xyz.y),  !  x, y, erg in m
!     *  sngl(aTrack.vec.w.x),  ! direc. cos.x in the current detector system.
!     *  sngl(aTrack.vec.w.y),  ! direc. cos.y
!     *  sngl(aTrack.vec.w.z),  ! direc. cos.z
!     * sngl(-angle.r(3)) ,
!     * sngl(aTrack.vec.coszenith),  ! cos of zenith angle
!     * sngl(inci.p.fm.p(4)-inci.p.mass)/Code2massN(int(inci.p.code)),
!     * inci.p.code
!         if(aTrack.p.code .eq. kelec) then
!            write(*, *) aTrack.where
!         endif
!      endif
!         you may need in some case other information such as
!     *  aTrack.p.subcode   ! sub code of the particle integer*2
!       aTrack.p.mass      ! mass 
!       aTrack.wgt         ! weight of the particle (may not be 1. if
!                           ! ThinSampling =T)
!       aTrack.p.fm.x      ! momentum x component.  Note. Momentum is
!                            given in the  Earth xyz system.

!       aTrack.p.fm.y      !          y
!       aTrack.p.fm.z      !          z
!        if(aTrack.p.code .eq. kelec .or. aTrack.p.code .eq. kphoton) 
!     *     then
!           ng = ng+1
!           sumg = sumg + aTrack.p.fm.e

!        endif
      endif
      end
      subroutine cwhere2dep(where, depth)
      implicit none
#include "Zcoord.h"
#include "Zpos.h"
#include "Zmagfield.h"
#include "Zobs.h"
#include "Zobsv.h"
 
      integer*2 where
      real*8 depth
      if(where .le. 0 .or. where .gt. NoOfSites) then
         write(*,*) where, NoOfSites
         stop 'error'
      endif

      depth =ObsSites(where).pos.depth
      end
!    *********************************** hook for end of 1 event
!    * At this moment, 1 event generation has been ended.
!    *
      subroutine chookEnEvent

      implicit none
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"




      type(track):: inci
      type(coord):: angle, tetafai
      integer i, j

!          ///////////
      integer nl, nth
      parameter (nl = 20, nth=12)
      common /testcos/ eth(nth),
     * ng(nth, nl), ne(nth, nl), nmu(nth, nl), 
     * nh(nth, nl), ntha
      integer ng, ne, nmu, nh, ntha
      real*8 eth
!     //////////////



!          ///////////
      real*8  fdepth, bsin, teta, fai, sumsize
      real*8  cgetBsin, sumx, sumy
      real*8 avex,  avey, sume
      integer nnew
!     //////////////
!       call cqIncident(inci, angle)
!      write(*,*) inci.vec.coszenith, angle.r(3)

      if(ObserveAS) then
         call cqFirstID(fdepth)
         fdepth = fdepth * 0.1         ! in g/cm2
         call cqIncident(inci, angle)
         angle.r(1) = -angle.r(1)   ! angle is directed to downward
         angle.r(2) = -angle.r(2)
         angle.r(3) = -angle.r(3)
         call ceCent2sph(angle, tetafai)
         teta = tetafai.r(1) 
         fai = tetafai.r(2)
         if(fai  .lt. 0. ) fai = 360.d0+fai
         bsin = cgetBsin(inci.p, Mag)*1.e4
!                   electron size in B approx.
!         write(*, *) (ASObsSites(i).esize, i=1, NoOfASSites)
!                   size weighted age
!         write(*, *) (ASObsSites(i).age,   i=1, NoOfASSites) 
         sumsize = 0.
!         write(*, *)
        do j = 1, ntha
          if(ntha .gt. 1) then
!             write(*,*) j
          endif
           do i = 1, NoOfASSites
              sumsize = sumsize + ASObsSites(i).esize
               write(*, '(f7.1,g13.3,f8.3,f7.1,
     *             4i8,f7.4)')
!     *       f8.3,    g13.3,f10.3,f10.3) ')
!     *           sngl(ASObsSites(i).pos.depth/10./angle.r(3)),
     *           sngl(ASObsSites(i).pos.depth/10.),
     *           sngl(ASObsSites(i).esize), 
     *           sngl(ASObsSites(i).age), sngl(fdepth),
!     *           sngl(bsin), sngl(sumsize), sngl(teta), sngl(fai)  
     *           ne(j, i), nmu(j, i), nh(j, i), ng(j, i), eth(j)
            enddo
         enddo
      endif

      end
!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun

      implicit none
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
!  #include  "Ztrackp.h"
      
      integer never   ! input & output
      
!         don't make never = 1, if you want to get
!         information after an electron made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the electron that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting
      never = 1
!
!        IntInfArray(ProcessNo).process will have one of
!       'brems', 'mscat', 'bscat', 'anihi' or 'mbrem'
!
      end

!     ********************* this is the hook called when
!       a gamma ray made an interaction.
!
      subroutine chookGInt(never)
            implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
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
!
!        default setting
      never = 1
!         IntInfArray(ProcessNo).process will have one of
!        'pair', 'comp', 'photoe' 'photop' 'mpair'
!       
      end

!     ********************* this is the hook called when
!       non e-g particle made an interaction.
!
      subroutine chookNEPInt(never)
            implicit none


#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"

!  #include  "Ztrackp.h"
      
      integer never   ! input & output
      
!         don't make never = 1, if you want to get
!         information after a non-e-g particle  made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the particle that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting
!
!      never = 1
        never = 1
        if(MovedTrack.p.code .eq.  kpion .or.
     *       MovedTrack.p.code .eq.  kkaon) then
            if(IntInfArray(ProcessNo).process .eq. 'coll') then
!              write(*,*)
!     *         MovedTrack.p.code,
!     *         sngl(MovedTrack.p.fm.p(4)), Nproduced
           endif
        endif
!
!        IntInfArray(ProcessNo).process  will have
!             'col' or 'decay'
      end

      





