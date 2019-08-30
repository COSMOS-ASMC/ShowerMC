#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
#if defined (KEKB) || defined (KEKA)
#define DOMPI
#endif
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      implicit none
#include "Zmanagerp.h"
 
!         If you feel writing the parameters on stderr is
!         a bother, comment out the next or
!         use other device than ErrorOut.
!         Also you may comment out all output routines below.

!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)


      end


!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      implicit none
#if defined (DOMPI)
#include "mpif.h"
#include "Zmpi.h"
      real*8 etime1
      etime1 = MPI_WTIME()
      write(0,*) 'rank=',mpirank, ' time start=',etime1
#endif
      end
  

!     ************************************ hook for observation
!     *  One particle information is brought here by the system.
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
#include "Ztrack.h"
      integer id  ! input.  1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
      if(id .eq. 2) then
!            output typical quantities.
        write(*, '(3i3, 1p, 4g15.4,0p, 4f10.6)')
     *  aTrack%where,   !  observation level. integer*2.  1 is highest.
     *  aTrack%p%code,    !  ptcl code.  integer*2.
     *  aTrack%p%charge,  !  charge,  integer*2 
     *  aTrack%t, !  relateive arrival time in nsec (NOT sec).
                        !  if TimeStructure is F, nonsense.
     *  aTrack%p%fm%p(4),  ! total energy in GeV.
     *  aTrack%pos%xyz%r(1), aTrack%pos%xyz%r(2), !  x, y in m
     *  aTrack%vec%w%r(1),  ! direc. cos%x in the current detector system.
     *  aTrack%vec%w%r(2),  ! direc. cos%y
     *  aTrack%vec%w%r(3),  ! direc. cos%z
     *  aTrack%vec%coszenith ! cos of zenith angle
      endif
!         you may need in some case other information such as
!       aTrack.p.subcode   ! sub code of the particle integer*2
!       aTrack.p.mass      ! mass 
!       aTrack.wgt         ! weight of the particle (may not be 1. if
                           ! ThinSampling =T)
!       aTrack.p.fm.p(1)      ! momentum x component.  Note. Momentum is
!                            given in the  Earth xyz system.

!       aTrack.p.fm.p(2)      !          y
!       aTrack.p.fm.p(3)      !          z

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

      integer i

      type(coord):: angle
      type(track):: incident
      
#if defined (DOMPI)
#include "mpif.h"
#include "Zmpi.h"
      real*8 etime2
      etime2 = MPI_WTIME()
      write(0,*) 'rank=',mpirank, ' time end=',etime2
#endif
      real(8):: xs
      integer:: IA, IZ 

      if(ObserveAS) then
!                   electron size in B approx.
         write(*, *) (ASObsSites(i)%esize, i=1, NoOfASSites)
!                   size weighted age
         write(*, *) (ASObsSites(i)%age,   i=1, NoOfASSites) 
      endif
      call cqIncident(incident, angle)

      write(0,*) '1ry c,subc,chg, TE   cos-zenith '
      write(0,*) 
     *   incident%p%code, incident%p%subcode, incident%p%charge,
     *     incident%p%fm%p(4), incident%vec%coszenith
      call cqFirstIPI(incident) ! first int point info.
      write(0,*) ' 1st int. depth, height, coszenith'
      write(0,*) incident%pos%depth, incident%pos%height,
     *     incident%vec%coszenith
      call cqFirstColMedia(IA, IZ, xs)
      write(0,*) ' first A,Z,xs=', IA, IZ, xs
      end


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun

      implicit none
#include "Ztrack.h"
#include "Ztrackp.h"
      integer klena
      character*24  tracefile
      character*1  qm/"'"/
      write(*,'(a)') '# ^             ^             ^'
      write(*,'(a)') '# |             |             | ' 
      write(*,'(a)') '# depth indx  code charge time   E '//
     *               ' X (m)     Y (m)       Wx       Wy  Wz'//
     *               ' cos(zenith)'

      write(*, *)

      write(0, *)
     * '        ****** Congratulations ****** ' 

      if(Trace .gt. 0) then
         tracefile = TraceDir(1:klena(TraceDir))//'trace?'
         write(0, *)
     * '       particle trace data has been created'//
     *     ' in '//tracefile
         write(0, *)
     * ' ?=1,2...   You can see it by Geomview (or gnuplot): '
         write(0, *)
     * 'Geomview is faster & better; If it has been installed, do '
         write(0, *)
     * '   dispcosTraceByGeomv  path-to-tracefile(=trace1 or 2)'
         write(0, *) "   "
         write(0, *)
     *     ' For gnuplot, use dispcosTraceByGnup instead '
         write(0,*) ' In both case, if no argument is given,'
         write(0,*) '  usage will be shown'
         write(0,*)
         write(0,*)
     * " ***Make Trace=0 in 'param file' for normal jobs***"  
      endif       
        write(0, *)
     * '************      Have a nice day !!      **************'

      end
!     ********************************* hook for trace
!     *  This is called only when trace > 60
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
!    if trace > 60. 
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

      h1 = TrackBefMove%pos%height- ObsSites(NoOfSites)%pos%height
      h2 = MovedTrack%pos%height - ObsSites(NoOfSites)%pos%height

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
!       'brems', 'mscat', 'bscat',or  'anihi'
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
!        'pair', 'comp', 'photoe' or 'photop'
!       
      end

!     ********************* this is the hook called when
!       non e-g particle made an interaction.
!

      subroutine chookNEPInt(never)
      use modXsecMedia
      implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
      
      integer never   ! input & output
      integer:: A,Z
      real(8):: xs, pabs, teta
      integer:: i



!         don't make never = 1, if you want to get
!         information after a non-e-g particle  made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the particle that made interaction
!   Pwork contains produced particles.
!   Nproduced has the 0number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting
      never = 1
      end

      
      

