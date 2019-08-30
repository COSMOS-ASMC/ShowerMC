#undef  USEHISTO
!    You can define USEHISTO only for intel fortran; 
!    (MacIFC, PCLinuxIFC)
!    You have to 'make' at UserHook/Hist for this
!    If make USEHISTO defined, interfaceHxxx.f's will be used for
!    taking  histogram on fly.
!
!=================================================================
      subroutine xBgRun
      implicit none
#include  "Zmaxdef.h"
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Ztrackp.h"
#include "Zcode.h"
#include "Zheavyp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include  "Zstackv.h"
#include  "Zincidentv.h"

#include "Zprivate.f"
#ifdef USEHISTO
#include "Zprivate1.f"
#include "Zprivate2.f"
#include "Zprivate3.f"     
      save rspec, lossrspec, arspec,  respec
      save rzspec,  zfspec, rtspec1, rtspec2, retspec1, retspec2
      save rezspec,  rzfspec, rfspec, efspec, refspec
#endif

      type(track):: incident
      type(coord):: AngleAtObs

      logical  HEobs            ! if T, currently observing 
      common /ZHEobs/ HEobs     !  particle is the one  obsrved at skeelton making time

            


      integer id  ! input.  1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack

      type(track):: inci
      type(coord):: angle
      type(coord):: tetafai
      
      character*128 input
      character*64 dirstr
      real sr, dr
      integer i, j, k, icon, ir, l
      integer NN
      integer iij, code, codex
      integer i1, i2, ic

      real*8 r, Eloss, rinmu, cosang
      real*8 dedt, rho, dist, disto, dedxmu
      real*8 fai, rminm

      real*8 u
      logical accept
      real*8 wx, wy, wz,  temp
      real   za
      real  de, Ek, f, molu
      real*8  cvh2den
      integer ldep, ridx, faiidx
      integer lengenv
      integer  ncpu ! # of smashed skeletons
      integer  mcpu ! and  skeletons to be used actully
      integer  margin ! # of additional skeletons fleshed for safety
      real*4 enhance  ! since we use only mcpu, the result must be enhanced 
                      ! by a factor of ncpu/mcpu
      integer binw
      character*9 ptcl(4)
      data ptcl/"Photons", "Electrons","Muons", "hadrons"/
      character*9 ptcl2(3)
      data ptcl2/"Electrons", "Muons","All"/
      real power(3)
      integer nstr
      real E0,  cosz, limit, limite
      data power/2.,2.,1./
      real  power2(3)
      data power2/2.,1.,2./
      character*128 title
      character*96 evid(nsites)
      real*8 cog, cog2, sumne,  obstimes, Savederg(5)
      real*8 firstz, dd
      logical dosort
      character*2  kd(3)
      integer kgetenv2
      integer leng,  lengn, lengid
      character*5 numb
      character*64 execid
      character*128 msg
                        
      save

#ifdef  USEHISTO
#include "interfaceHBGR.f"
#endif
!             used only if enhace >1; but the value may be different
      limit = 20000
      lengenv = kgetenv2("NCPU", input) 
      read(input(1:lengenv),*)  ncpu
      lengenv = kgetenv2("MCPU", input)
      read(input(1:lengenv),*)  mcpu
      lengenv = kgetenv2("MARGIN", input)
      read(input(1:lengenv),*)  margin
      enhance = ncpu
      enhance = enhance/mcpu
      limite = limit * enhance  ! this must be enhanced from the first.
      write(0,*) ' *** ncpu=',ncpu,
     *            ' mcpu=',mcpu, ' enhance factor=',enhance
      write(0,*) ' margin=',margin
      input = ' '
      lengn =  kgetenv2("NUMB", numb)
      leng =  kgetenv2("OUTDIR", input)
      lengid = kgetenv2("EXECID", execid)
      
      if(ObserveAS) then
         msg =input(1:leng)//"/"//execid(1:lengid)//
     *        "-@."//numb(1:lengn)//".hyb"
         call copenfw2(fnoB, msg, 1, icon)
         if(icon .gt. 1) then
            write(0,*) ' icon=', icon
            call cerrorMsg(msg, 1)
            call cerrorMsg('could not be opened', 0)
         endif
      endif
      return
!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      entry xBgEvent

      obstimes = 0.
      call cqIncident(inci, angle)

      cosz = -angle.r(3)
      E0 = inci.p.fm.p(4)
      if(inci.p.code .eq. 9) then
         NN= inci.p.subcode
      elseif(inci.p.code .eq. 1) then
         NN=0
      else
         NN=1
      endif
!         next is only available for parallel job.  For normal job,
!         fisrt col.depth is not yet fixed.
      firstz= Zfirst.pos.depth*0.1
      write(0,'(a,1pE11.3,a,E11.3,a,E11.3,a)')
     *     ' E0=',E0, ' cosz=',cosz, ' firstz=',
     *    firstz, ' g/cm2' 
!      
      do i = 1, NoOfASSites
         SumEloss(i) = 0.
         Ng(i) = 0.
         Ne(i) = 0.
         Nmu(i) = 0.
         Nhad(i) = 0.
      enddo
#ifdef USEHISTO
#include "interfaceHBGE.f"
#endif
      return
!     ***************
      entry xObs(aTrack, id)
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
!     **************************
!            to be able to see the job is really running (i.e, no loop)
!         we write messages every 500000 calls.
!
      obstimes = obstimes + 1.d0
      if(mod(obstimes, 500000.d0) .eq. 0. ) then 
         dosort=.false.
         do i = 1, min(4,Stack_pos)
            if(Stack(i).p.fm.p(4) .ne. Savederg(i)) then
               Savederg(i)=Stack(i).p.fm.p(4) 
               dosort=.true.
            endif
         enddo
         if(dosort) then
            call csortStack
         endif
         write(0, *) ' obstimes=', obstimes, ' ptclE=',aTrack.p.fm.p(4)
         do i = 1, min(4,Stack_pos)
            write(0,*)' stack tops=', Stack(i).p.fm.p(4)
         enddo
      endif
!     ***************
      code = aTrack.p.code
      if(id .eq. 2 .and. code .le. 6 ) then
         codex=min(code, 4)
         wz = aTrack.vec.w.r(3) ! downgoing < 0
         if(wz .gt. 0) return
         wz = -wz
         ldep =  aTrack.where

         if(code .eq. kphoton) then
            Ng(ldep) = Ng(ldep) + enhance
         elseif(code .eq. kelec) then
            Ne(ldep) = Ne(ldep) + enhance
         elseif(code .eq. kmuon) then
            Nmu(ldep) = Nmu(ldep) + enhance
         else
            Nhad(ldep)=Nhad(ldep) + enhance
         endif

         Ek = aTrack.p.fm.p(4) -aTrack.p.mass

!            ---------- compute energy loss rate
         if(aTrack.p.charge .ne. 0  ) then
            rho = cvh2den(aTrack.pos.height)
!         get energy loss when aTrack goes 1 g/cm2 along the
!         primary direction. Gramage the particle can run is
!         1/cos where cos is the cos of angle relative to the
!         primary angle . 1g/cm^2 = 10-3kg/10-4 m^2 =10 kg/m^2.
!         To travel  1 g/cm^2  along shower axis, the ptcl must
!         run dist kg/m^2
            if(abs(wz) .gt. 1.d-2) then
               dist =10./wz     ! in kg/m2/(g/cm2)
            else
!                      for safety
               dist =1000.
            endif
            if( HEobs ) then
!                   the ptcls is the one obsrved at skeleton making time
!                   we must compute dedt here
               call cdedxInAir(aTrack.p, rho, dedt) ! dedt; GeV/(kg/m2)
               if(aTrack.p.code .eq. kmuon ) then
!                dE/dx due to muon pair, brem, nuc.i
                  call cmudEdx(MuNI, MuBr, MuPr, aTrack.p.fm.p(4),
     *                 dedxmu)  ! dedxmu in GeV/(g/cm2)
                  dedxmu = dedxmu /10. !  GeV/(kg/m2)
                  dedt = dedt + dedxmu
               endif
            else
!               we can use already computed one
               call cqElossRate(dedt) !  loss rate GeV/(kg/m^2)
            endif
!                       energy loss rate
            Eloss = dedt*dist*enhance !  GeV/(g/cm2)

            SumEloss(ldep)=SumEloss(ldep) + Eloss
         else
            Eloss=0.
         endif
         if(aTrack.where .eq. 6) then
            write(*,
     *           '(4i3, 1p4E11.3, 0p, 2f8.4,f10.6)')
     *        ldep,  code,  aTrack.p.subcode, aTrack.p.charge, 
     *        Ek, aTrack.t, 
     *        aTrack.pos.xyz.x, aTrack.pos.xyz.y,
     *        -aTrack.vec.w.r(1),  -aTrack.vec.w.r(2),  wz
         endif
      endif
#ifdef USEHISTO
#include "interfaceHOBS.f"
#endif
      return
!     **************
      entry xEnEvent
!     **************
      firstz= Zfirst.pos.depth*0.1
      if(ObserveAS) then
         cog = 0.
         sumne = 0.

         do i = 1, NoOfASSites
            ASObsSites(i).esize = ASObsSites(i).esize* enhance
            if(i .gt. 1 .and. i  .lt. NoOfASSites ) then
               dd =(ASDepthList(i+1) - ASDepthList(i-1))/2.0
            elseif(i .eq. 1) then
               dd =(ASDepthList(2) - ASDepthList(1))
            else
               dd =(ASDepthList(NoOfASSites) -
     *              ASDepthList(NoOfASSites-1))
            endif
            cog = cog + ASObsSites(i).esize*dd*ASDepthList(i)
            sumne= sumne +ASObsSites(i).esize*dd
         enddo
!          0.1 is for g/cm2
         cog = cog*0.1/sumne
         cog2 = 0.
         sumne = 0.
         do i = 1, NoOfASSites
            if( ASObsSites(i).age .gt.
     *          (2.0-ASObsSites(NoOfASSites).age))  then
               if(i .gt. 1 .and. i  .lt. NoOfASSites ) then
                  dd =( ASDepthList(i+1) - ASDepthList(i-1))/2.0
               elseif(i .eq. 1) then
                  dd =(ASDepthList(2) - ASDepthList(1))
               else
                  dd =(ASDepthList(NoOfASSites) -
     *              ASDepthList(NoOfASSites-1))
               endif
               dd = dd
               cog2 = cog2 + ASObsSites(i).esize*ASDepthList(i)*dd
               sumne= sumne +ASObsSites(i).esize*dd
            endif
         enddo
         if(sumne .gt. 0.) then
            cog2 = cog2*0.1/sumne
         else
!              too deep penetration
            cog2 = ASDepthList(NoOfASSites)*0.1
         endif

         if(fnoB .ge. 0 )  then
            write(fnoB,
     *      '("h ", i4,  3i3, 1pE11.3, 0p 3f11.7, f7.2, 2f7.0)')
     *      EventNo,  inci.p.code,
     *      inci.p.subcode, inci.p.charge,
     *      inci.p.fm.e, -angle.r(1), -angle.r(2), -angle.r(3),
     *      firstz, cog, cog2
         endif
         do i = 1, NoOfASSites 
            if(fnoB .ge. 0) then
               write(fnoB, '("t ", i3, 2f7.1,  2f6.3,
     *         1p6E11.3)')
     *          i, 
     *          ASDepthList(i)*0.1,  ASObsSites(i).mu,
     *          ASObsSites(i).age,   ASDepthList(i)*0.1/cog2, 
     *          Ng(i), Ne(i), Nmu(i), Nhad(i),
     *          ASObsSites(i).esize, SumEloss(i)  
            endif
         enddo
         if(fnoB .gt. 0 ) then
            write(fnoB,*)
         endif
      endif
#ifdef USEHISTO
#include "interfaceHENE.f"
#endif
      end
