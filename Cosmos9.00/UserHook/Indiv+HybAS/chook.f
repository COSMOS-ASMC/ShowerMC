!     take energy spectrum of e, g, mu, had, neumu, neue
!     at different depths   (going down and going up)      
!     parametarizing  by r
#include "cmain.f"
#include "chookHybAS.f"
#include "ctemplCeren.f"
      module modASsize
      use modhistogram2
      implicit none
      integer,save::Nlayers=0    ! no of depths
      integer,save::nev=0       ! event counter
      integer,parameter:: nbinR=2*4 ! one log10 -->2.  1 to 10^4 m

! take energy spectrum with R being a parameter
!     at the deepest layer or a layer specified by the next (if>0)
      integer,save:: nthLayer=0
      integer,parameter:: nkind=6 ! particle type classification
                !  g, e, mu, (pi,K,N), neue, neumu
      integer,parameter:: nbinE=20*nkind ! one log10 -->20.  6 decades
                                !    for g,e. >=100 keV  K.E
                                !    for mu,had, neu >=10 MeV K.E
      real(4),save::binE=1./(nbinE/nkind)
      real(4),save::binR=1./(nbinR/4)
      !                   g,e,mu,had,neue,neumu
      type(histogram2):: ERhist(nkind), TERhist(nkind)

      real(4),save::EminHist(nkind)
      data EminHist /2*100e-6,4*30e-3/
      
!     transition curve
      integer,parameter::nEmin=2
      integer,parameter::nkindTran=7
      !          1 2   3  4            5    6      7 
      !          g,e, mu,(pi,K,N,A), neue, neumu, others
      real(8),save::EminTran(nkindTran, nEmin)
      data EminTran/2*100.d-6,5*30d-3, 2*500d-6, 5*1.d0/
!     Nlayer,nkindTran, down-or-up
!     g,e,mu, (had), neue, neumu, others==>7
      ! first is layer, 2nd emin, 3rd kind, 4th updown
      real(8),allocatable,save::Num(: ,:,:,:)
      ! first layer, 2nd emin, 3rd kind
      real(8),allocatable,save::pnNum(:,:,:)  ! p,n only down going
      integer,parameter::fno=21
      end module modASsize
      
!  *************************************** hook for Beginning of a Run
!  * At this moment, all (system-level) initialization for this run
!  * has been ended.  After this routine is executed, the system goes into the
!  * event creation loop.
!  *
      subroutine chookBgRun
      use modASsize
      use modHistogram1
      implicit none
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"

!         If you feel writing the parameters on stderr is
!         a bother, comment out the next or
!         use other device than ErrorOut.
!     Also you may comment out all output routines below.

      integer:: i
      
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
      if(NoOfSites /= NoOfASSites) then
          write(0,*)
     *    '# of Depth and ASdepth in this appli. must be the same'
          stop
      endif
      Nlayers = NoOfSites
      if( any( ASObsSites(1:Nlayers).pos.depth /=
     *         ObsSites(1:Nlayers).pos.depth ) ) then
         write(0,*)
     *  'Depth and ASDepth must be the same in this appli.'
         stop
      endif
      if( nthlayer > Nlayers ) then
         write(0,*)' nthlayer > # of given depths=',Nlayers
         stop
      elseif( nthlayer == 0 ) then
         nthlayer = Nlayers
      endif
         
         
      open(fno, file='try.hist', form='formatted')
      call kwhistso(1)   ! ascii output 
      do i =1,  nkind
         !  instance of each ptcl class: kind
         call kwhisti2(ERhist(i),  
     *     EminHist(i), binE, nbinE, b'00011',
     *        1.0,  binR, nbinR, b'00001' )
         call kwhisti2(TERhist(i),  
     *     EminHist(i), binE, nbinE, b'00011',
     *        1.0,  binR, nbinR, b'00001' )
      enddo
      !  clear the histogram area for total events
      do i=1, nkind
         call kwhistc2(TERhist(i))
      enddo
      ! counter for transition of nkindTran 
      allocate( Num(Nlayers, nEmin, nkindTran, 2) )

!     2 is for p,n (pbar,nbar included); specially take
      !   p and n tran. only down goingn
      allocate(pnNum(Nlayers, nEmin, 2) )
      end
#ifdef sun4
      integer function csigHandler(sig, code, context)
      implicit none
#include "Zmanagerp.h"
      integer sig, code, context(5)
      write(ErrorOut, *)  ' f.p exception content=' , context(4)
!      call abort()
      end function csigHandler
#endif

!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      subroutine chookBgEvent
      use modASsize
      implicit none
!          ///////////
#include "Ztrack.h"
#include "Zobs.h"
#include "Zobsv.h"


      integer i

      integer seed(2)
!      write(*, *) ' bigin event generation'
      call cqIniRn(seed)
!     write(*,*) ' seed=', seed
       !  clear the histogram area of this event
      do i=1, nkind
         call kwhistc2(ERhist(i))
      enddo

!          clear the transtion counter of this event
      Num(:,:,:,:) = 0.
      pnNum(:,:,:) = 0.
      end subroutine chookBgEvent
  

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
      use modASsize
      implicit none
#include "Zcode.h"
#include "Ztrack.h"
      integer id  ! input.  1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track)::aTrack

      integer:: kind, updw, iemin, kindTran
      real(4)::KE,R
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
      if(id .eq. 2) then
!            output typical quantities.
!        write(*, *) 
!     *  aTrack.where,   !  observation level. integer*2.  1 is highest.
!     *  aTrack.p.code,    !  ptcl code.  integer*2.
!     *  aTrack.p.charge,  !  charge,  integer*2 
!     *  sngl(aTrack.t), !  relateive arrival time in nsec (NOT sec).
!c                        !  if TimeStructure is F, nonsense.
!     *  sngl(aTrack.p.fm.e),  ! total energy in GeV.
!     *  sngl(aTrack.pos.xyz.x), sngl(aTrack.pos.xyz.y), !  x, y in m
!c     *  sngl(aTrack.vec.w.x),  ! direc. cos.x in the current detector system.
!c     *  sngl(aTrack.vec.w.y),  ! direc. cos.y
!c     *  sngl(aTrack.vec.w.z),  ! direc. cos.z
!     *  sngl(aTrack.vec.coszenith) ! cos of zenith angle
!c         if(aTrack.p.code .eq. kelec) then
!            write(*, *) aTrack.where
!         endif
!      endif
!         you may need in some case other information such as
!       aTrack.p.subcode   ! sub code of the particle integer*2
!       aTrack.p.mass      ! mass 
!       aTrack.wgt         ! weight of the particle (may not be 1. if
!                           ! ThinSampling = T)
!       aTrack.p.fm.x      ! momentum x component.  Note. Momentum is
!                            given in the  Earth xyz system.

!       aTrack.p.fm.y      !          y
!       aTrack.p.fm.z      !          z
         kind=aTrack.p.code
         
         KE = aTrack.p.fm.p(4) - aTrack.p.mass
         do iemin = nEmin, 1, -1
            if( KE >= EminTran(min(kind,9),iemin) ) goto 10
         enddo
         iemin= 0
 10      continue
         if( iemin > 0  ) then
            if( kind == 6 .and. aTrack.vec.coszenith >0.) then
               ! p,n (or pbar,nbar)
               if( aTrack.p.charge == 0 ) then
                  pnNum(aTrack.where, iemin, 2) =
     *                 pnNum(aTrack.where, iemin, 2) + aTrack.wgt
               else
                  pnNum(aTrack.where, iemin, 1) =
     *                 pnNum(aTrack.where, iemin, 1) + aTrack.wgt
               endif
            endif
            if( kind <= 3 ) then
               kindTran= kind
            elseif(kind >= 4 .and. kind <=6 .or. kind == 9 ) then
               kindTran = 4         !  pi,K, n,A
            elseif( kind == 7 ) then
               kindTran = 5  ! 7,8  neue
            elseif( kind ==8 ) then
               kindTran = 6    ! 8  neumu
            elseif(kind > 8) then
               kindTran = 7     ! others
            endif
            if(aTrack.vec.coszenith > 0.) then
               updw = 1
            else
               updw = 2
            endif

            Num(aTrack.where, iemin, kindTran,  updw) =
     *       Num(aTrack.where, iemin, kindTran, updw) + aTrack.wgt
         endif
!-------
         if( aTrack.where == nthlayer ) then
            if(aTrack.vec.coszenith > 0.) then
               if(kindTran < 7 ) then
                 R = sqrt(aTrack.pos.xyz.x**2 + aTrack.pos.xyz.y**2)
                 call kwhist2( ERhist(kindTran),  KE, R,
     *                  sngl(aTrack.wgt) )
               endif
            endif
         endif
      endif
      end

!    *********************************** hook for end of 1 event
!    * At this moment, 1 event generation has been ended.
!    *
      subroutine chookEnEvent
      use modASsize
      implicit none
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"




      integer:: i, j

      nev= nev + 1
      do i = 1, nkind 
         if( nev == 1 ) then
            TERhist(i) = ERhist(i)
         else
            call kwhista2(TERhist(i), ERhist(i), TERhist(i)) ! accumulation
         endif
      enddo

!----------transition: Num; Emin matter
      do i = 1, Nlayers
         do j = nEmin-1, 1, -1
            Num(i, j, :, :) = Num(i, j, :, :) + Num(i,j+1,:,:)
            pnNum(i, j, :) =pnNum(i,j,1) + pnNum(i, j+1, :)
         enddo
      enddo
      
      if(ObserveAS) then
!     electron size in B approx and size weighted age
          ! counter for transition of nkindTran 
          !   Num(Nlayers, nEmin, nkindTran, 2) 
         write(*,'(a)')
     *   "# id L  v dep     H    age   Nehyb        Ne1        Ne2"
         do i = 1, Nlayers
            write(*, '(a, i2, f8.1, f8.0, f5.2, 1p, 3g12.5)')
     *            "Ne ", i, ObsSites(i).pos.depth/10.,
     *      ObsSites(i).pos.height, ASObsSites(i).age,
     *      ASObsSites(i).esize, Num(i, 1, 2, 1), Num(i, 2, 2, 1)
         enddo
      endif
!-----------n,p
      write(*,'(a)')
     *     "# id  L   v dep    H      Np1         Nn1 "//
     *     "      Np2         Nn2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        "Np/n ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height, pnNum(i, 1, 1), pnNum(i,1, 2),
     *        pnNum(i, 2, 2), pnNum(i, 2, 2)
      enddo
!-------------  g tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Ngd1        Ngu1 "//
     *     "       Ngd2       Ngu2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        " Ng ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 1, 1), Num(i, 1, 1, 2),
     *        Num(i, 2, 1, 1), Num(i, 2, 1, 2)              
      enddo

!-------------  e tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Ned1        Neu1 "//
     *     "       Ned2        Neu2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        " Ne ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 2, 1), Num(i, 1, 2, 2),
     *        Num(i, 2, 2, 1), Num(i, 2, 2, 2)              
      enddo

      !-------------  mu tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Nmd1        Nmu1 "//
     *     "       Nmd2        Nmu2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        "Nmu ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 3, 1), Num(i, 1, 3, 2),
     *        Num(i, 2, 3, 1), Num(i, 2, 3, 2)              
      enddo

!-------------  had tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Nhd1        Nhu1 "//
     *     "       Nhd2        Nhu2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        " Nh ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 4, 1), Num(i, 1, 4, 2),
     *        Num(i, 2, 4, 1), Num(i, 2, 4, 2)              
      enddo

      !-------------  nue tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Nned1       Nneu1 "//
     *     "       Nned2        Nneu2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        "Nnue ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 5, 1), Num(i, 1, 5, 2),
     *        Num(i, 2, 5, 1), Num(i, 2, 5, 2)              
      enddo


!-------------  numu tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Nnmd1       Nnmu1 "//
     *     "       Nnmd2        Nnmu2"
      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        "Nnum ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 6, 1), Num(i, 1, 6, 2),
     *        Num(i, 2, 6, 1), Num(i, 2, 6, 2)              
      enddo


      !-------------  numu tran
      write(*,'(a)' )
     *     "# id L  v dep     H      Nod1       Nou1 "//
     *     "       Nod2        Nou2"

      do i = 1, Nlayers
         write(*,'(a, i2, f8.1, f8.0,  1p, 4g12.5)')
     *        "Noth ", i, ObsSites(i).pos.depth/10.,
     *        ObsSites(i).pos.height,
     *        Num(i, 1, 7, 1), Num(i, 1, 7, 2),
     *        Num(i, 2, 7, 1), Num(i, 2, 7, 2)              
      enddo

!       separator
      write(*, *)

      end      subroutine chookEnEvent


!     ********************************* hook for end of a run
!     *  all events have been created or time lacks
!     *
      subroutine chookEnRun
      use modASsize
      implicit none
      end       subroutine chookEnRun
!     ********************************* hook for trace
!     *  This is called only when trace > 100
!     *  User should manage the trace information here.
!     *  If you use this, you may need some output .for trace
!     *  at the beginning of 1 event generatio and at the end of  1 event
!     *  generation so that you can identfy each event.
!     *
!     *
      subroutine chookTrace
      use modASsize
      implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsv.h"
#include  "Zcode.h"
       type(coord)::f
       type(coord)::t
       logical compress/.true./
       save compress


!       real*4 h1,  h2
!
!    If trace > 100.
!    every time a particle is moved in the atmosphere, this routine is called.
!
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

!      h1 = TrackBefMove.pos.height- ObsSites(NoOfSites).pos.height
!      h2 = MovedTrack.pos.height - ObsSites(NoOfSites).pos.height
!      This example here is to put only muons with the same format
!       the standard trace infomation.  Trace-100 is treated as
!       the standard Trace value.
!     'compress' is to indicate compressed output.
!
!      if(TrackBefMove.p.code .eq. kmuon) then
!         call ccoordForTr(Trace-100, f, t)   ! convert coordinate
!         call cwrtTrInfo(compress, f, t)    ! write trace data
!      endif

      end      subroutine chookTrace

!     ********************* this is the hook called when
!       an electron made an interaction.
!
      subroutine chookEInt(never)
      use modASsize
      implicit none

#include  "Ztrack.h"
#include  "Ztrackv.h"
!  #include  "Ztrackp.h"
      
      integer never   ! input & output
      
!         don't make never = 1, if you want to get
!         information after an electr.on made interaction
!         if this is made non zero, this routine will never be called.
!
!   MovedTrack is the electron that made interaction
!   Pwork contains produced particles.
!   Nproduced has the number of particles in Pwork
!   IntInfArray(ProcessNo) contains the type of interaction
!
!        default setting
      never = 1
!       never = 2: same as 0
!             = 3: discard all the child of this (but not other ptcls)
!             = 4: discard the event
!
!        IntInfArray(ProcessNo).process will have one of
!       'brems', 'mscat', 'bscat', 'anihi' or 'mbrem'
!
      end      subroutine chookEInt

!     ********************* this is the hook called when
!       a gamma ray made an interaction.
!
      subroutine chookGInt(never)
      use modASsize
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
!       never = 2: same as 0
!             = 3: discard all the child of this (but not other ptcls)
!             = 4: discard the event

!         IntInfArray(ProcessNo).process will have one of
!        'pair', 'comp', 'photoe' ,'photop' or 'mpair'

!       
      end       subroutine chookGInt

!     ********************* this is the hook called when
!       non e-g particle made an interaction.
!
      subroutine chookNEPInt(never)
      use modASsize
      implicit none

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
      never = 1
!       never = 2: if proton is the current particle, give 2 
!                  if you make a hybrid air shower from that and
!                  want to discard that proton.
!             = 3: discard all the child of this (but not other ptcls)
!             = 4: discard the event

!
!        IntInfArray(ProcessNo).process  will have
!             'col' or 'decay'
      end subroutine chookNEPInt
      
