#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
!     If you would supply your own cmyEfield.f and/or cmyBfield.f
!     put your cmyEfield.f in this folder; to do so,
!     probably it's better to copy cmyEfield.f in $COSMOSTOP/cosmos/
!     here and modify it.  The file here will override the one
!     in $COSMOSTOP/cosmos/
#include "cmyEfield.f"
#include "cmyBfield.f"

#if defined (KEKB) || defined (KEKA)
#define DOMPI
#endif

      module mycounter
      implicit none
      integer,parameter:: nEmin=5
      real(8),parameter:: EminV(nEmin)= (/100.d-6, 300.d-6, 1.d-3,
     *     3.0d-3,10.d-3/)

                     
      real(8),allocatable:: Num(:, :, :)
      real(8),allocatable:: NumT(:, :, :)
      real(8),allocatable:: ASsize(:) ! where->i
      real(8),allocatable:: ASage(:) ! where->i
!          ==>   Num(nEmin, 6, NoOFSites)
      ! 6=(g,e,mu,pi,K,N )
      end module mycounter
      
      
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      use mycounter
      implicit none
#include "Zmanagerp.h"

#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

      real(8):: oldv
      integer:: icon
      
!         If you feel writing the parameters on stderr is
!         a bother, comment out the next or
!         use other device than ErrorOut.
!         Also you may comment out all output routines below.

!
!            namelist output
      call cwriteParam(ErrorOut, 0)
!            primary information
      call cprintPrim(ErrorOut)
!            observation level information
      call cprintObs(ErrorOut)
!     reset critical energy; this is used only for "Air"  and
!     "Generate" parameter contains "as". The deafult value (~86MeV) is
!     calculated as defined in PDB but for the Air shower size
!     caclucalation, it seems to large.
      call epResetEcrit(0, "Air", 81.0d-3, oldv, icon)
      write(0,*) 'icon=',icon,'Default Ecrit oldv(MeV)=', oldv*1000,
     * ' has been reset to ', 81, ' MeV'

      allocate( Num(nEmin,6, NoOfSites) )
      allocate( NumT(nEmin,6, NoOfSites) )      

      allocate( ASsize(NoOfASSites) )
      allocate( ASage(NoOfASSites) )
      NumT(:,:,:) = 0.d0
      end
#if defined DontDEFINE
!     This is to show examples of how to use the generic particle code
!     conversion routine:  cpdg2cos and ccos2pdg
!==================      
!     use modCodeConv  ! must be put before implicit none
!     #include "Zptcl.h"  ! if use type(ptcl) without type(track)
!      
!       type(ptcl):: aP, bP
!       integer:: code, subcode, charge
!       integer(2):: hcode, hsubc, hchg
!       integer:: pdgc, pdgc2
!
!          assmue nP is ready      
!        call ccos2pdg(aP, pdgc)
!          assume code, subcode, charge are ready
!        call cos2pdg(code, subcode, charge, pdgc)
!          pbar
!        call ccos2pdg(6, 1, -1, pdgc)
!          assume pdgc is ready ; get integer      
!        call cpdg2cos(pdgc, code, subcode, charge)
!         assume    pdgc is ready, get integer(2)
!        call cpdg2cos(pdgc, hcode, hsubc, hchg)
!         assume pdgc is ready, get in bP
!        call cpdg2cos(pdgc,bP)
!
#endif      
!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      use mycounter      
      implicit none
#include "Ztrack.h"
#include "Ztrackv.h"            
#include "Zobs.h"      
#include "Zobsv.h"      
      integer:: enum, cumnum
#if defined (DOMPI)
#include "mpif.h"
#include "Zmpi.h"
      real*8 etime1
      etime1 = MPI_WTIME()
      write(0,*) 'rank=',mpirank, ' time start=',etime1
#endif
      call cqEventNo(enum, cumnum)
      write(*,*) '## ev', enum
      Num(:,:,:) = 0.
      
      end
  

!     ************************************ hook for observation
!     *  One particel information is brought here by the system.
!     *  All information of the particle is in aTrack
!     *
      subroutine chookObs(aTrack, id)
      use mycounter
      use modCodeConv
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

      integer::i, code
      real(8):: KE

      integer::pdgcode
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
      if(id .eq. 2) then
         call ccos2pdg(aTrack%p, pdgcode)
         KE = aTrack%p%fm%p(4)-aTrack%p%mass
!            output typical quantities.
         write(*, '(4i4,i12, 1p, 4g15.4,0p, 4f10.6)')         
     *  aTrack%where,   !  observation level. integer*2.  1 is highest.
     *  aTrack%p%code,    !  ptcl code  integer*2.
     *  aTrack%p%subcode,    !  ptcl code  integer*2.
     *  aTrack%p%charge,  !  charge,  integer*2
     *  pdgcode,          ! kf (pdg) code
     *  KE,    ! Kinetic enrgy    
     *  aTrack%t, !  relateive arrival time in nsec (NOT sec).
                 !  if TimeStructure is F, nonsense.

     *  aTrack%pos%xyz%r(1:2), ! x,y in m
     *  aTrack%vec%w%r(1:3),  ! direc. cos.x,y,z in the current det. system.
     *  aTrack%vec%coszenith ! cos of zenith angle
         code=aTrack%p%code 
         if( code > 0 .and. code < 7) then
            do i = 1, nEmin
               if( KE >= EminV(i) ) then
                 Num(i, code, aTrack%where) =
     *                 Num(i, code,aTrack%where) +
     *                 aTrack%wgt
               else
                  exit
               endif
            enddo
         endif         
      endif

!        To convert the momentum into the observational
!       coordinate system,  you may call
!       call cresetMom(aTrack)
      end

!    *********************************** hook for end of 1 vent
!    * At this moment, 1 event generation has been ended.
!    *
      subroutine chookEnEvent
      use mycounter
      implicit none
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
      

      type(coord):: angle
      type(track):: incident
      real(8):: xs
      integer:: IA, IZ
      integer:: i, code
!      real(8):: h0, cos0, h, t0, len, tslan, tp, Rcos, dh, cterm
      real(8),external::  clen2thick, clen2thickEx !,cgetbsin,
      

#if defined (DOMPI)
#include "mpif.h"
#include "Zmpi.h"
      real*8 etime2
      etime2 = MPI_WTIME()
      write(0,*) 'rank=',mpirank, ' time end=',etime2
#endif
      
      call cqIncident(incident, angle)
      write(*,*) '#1ry c,subc,chg, TE   cos-zenith '
      write(*,'(3i4,1p, g14.4,0p, f9.5)')
     *     incident%p%code, incident%p%subcode, incident%p%charge,
     *     incident%p%fm%p(4), incident%vec%coszenith
      call cqFirstIPI(incident) ! first int point info.
      write(*,*) '#1st int. depth, height, coszenith'
        
      write(*,'(1p,3g14.4)')
     *     incident%pos%depth, incident%pos%height,
     *     incident%vec%coszenith
      call cqFirstColMedia(IA, IZ, xs)
      write(*,*) ' first col. target  A,Z,xs=', IA, IZ, xs


      if(ObserveAS) then

         ASsize(:) = ASsize(:) +     ASObsSites(i)%esize
         ASage(:) = ASage(:) +    ASObsSites(i)%age 

         do  i = 1, NoOfASSites         
            write(*,
     *           '(a, i3, f8.1, 1p, g14.4, 0p, f7.3)')
     *           '#HybAS  ',
     *           i, ASObsSites(i)%pos%depth/10.,
     *           ASObsSites(i)%esize,
     *           ASObsSites(i)%age
         enddo
      endif      

      do i =1, noOfSites
         do code =1,6
            write(*,
     *       '("#avNe (>thres)",i3, i3, 1p,10g14.4)')
     *           i, code, Num(1:nEmin, code, i)
         enddo
         NumT(:, code, i) = NumT(:,code, i) +
     *       Num(:,code, i)         
      enddo
         
         
!        ************ if you want to flesh this event later
!        you may keep the random no. seed  by the following
!        call cwriteSeed   !  SeedFile
      end


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun
      use mycounter
      implicit none
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
      
      

      integer:: enum, cumnum
      integer:: i, code

      call  cprintStatus  ! if don't like,  comment o
      call  cqEventNo(enum, cumnum)
!      write(*, '("#KE Emin(GeV)=", 1p,10g14.4)') EminV(1:nEmin)

      do i =1, NoOfSites
         do code =1,6
            write(*,
     *       '("#avN (>thres)",i, 1p,10g14.4)')
     *        i, code, NumT(1:nEmin, code, i)/enum
         enddo
      enddo
      

      if(ObserveAS) then      
         do  i = 1, NoOfASSites
            write(*,
     *         '("#AS size", i3, 1p, g14.4, 0p, f6.3)')
     *       i, ASsize(i)/enum, ASage(i)/enum
         enddo
      endif
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

      h1 = TrackBefMove%pos%height- ObsSites(NoOfSites)%pos%height
      h2 = MovedTrack%pos%height - ObsSites(NoOfSites)%pos%height

      end
!     ********************* this is the hook called when
!       an electron made an interaction.
!
      subroutine chookEInt(never)
            implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zpwork.h"            
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
!       'brems', 'mscat', 'bscat'  'anihi' or 'mbrem'
!
      end

!     ********************* this is the hook called when
!       a gamma ray made an interaction.
!
      subroutine chookGInt(never)
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zpwork.h"   

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
!         IntInfArray(ProcessNo)%process will have one of
!     'pair'  pair creation
!     'comp'  compton effect
!     'coh'   coherent scatterign  
!     'phot'  photo-eelcric effect  (old one was photoe)
!     'mpair' magnetic pair creation
!     'photop' photo-hadron production
!      
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
#include  "Zpwork.h"   
      
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
      never = 1
!
!     IntInfArray(ProcessNo)%process  will have
!     "knoc"  knock-on by charged partile (
!     "coll" (hadron collision)
!     "decay"  decay
!     "pair"   pair creation by muon
!     "brem"   brems by muon
!     "nuci"   hadron prodution by muon
!     "capt"   mu- capture
      end

      
