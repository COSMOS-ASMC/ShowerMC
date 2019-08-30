!      make next <=0 for ascii write for individual ptcl info.
!      Normmally binary output is gathered and converted into
!      ascii in the last stage so this is kept as it is.
!      33 is used as fortran file no.
#define FNODATDEF 33
!          standard
#define BUFSIZE  100000
!          prog. size upto 512MB at kekb.  no disk I/O
!  #define BUFSIZE 2000000
!       define next if reduce the output size(.dat and rec of .nrfai)
!       on fly even at non-mpi env.

!
#if defined (KEKB) || defined (KEKA) 
#define DOMPI
#endif
!   output:
!    Layer specification where output data is taken is
!    determined as follows
!    A) ASDepthList.   This is for .hyb data.  
!         ASDepthList(i), i=1, NoOfASSites
!    B) DepthList.  This is for web "all" and web "dEdx"
!         DepthList(i), i=1, NoOfSites
!        B') With additional param.   
!           indivdep(i)  for individual ptcl output.
!           (i=1, ansites) (ansites<=NoOfSies). j=indivdep(i)
!           gives the layer number.
!           histdep(i)   for histogram  output.
!           (i=1, hnsites) (hnsites<=NoOfSites). j=histdep(i)
!           gives the layer number.
! 
!    1)   web data
!            rec: count no. of particles of each type
!                 for which we recorded individual ptcl
!                 type info.   (indivdep)
!         
!            all: for all layers. count the number of
!                     each ptcl type which passes the
!                     given layer.  
!           dEdx: for all layers. count dE/dx of charged
!                 ptcls which pass the given layer
!                  
!    2)  hyb data
!               all layers (NoOfASSites)
!    3) .dat   for only specified 1 or few layers. individual
!              info.  corresponding to web "rec" data.
!    4) .hist  sub set of .dat layers.
!          lateral distribution corresponding to each fai
!             bin of the web sector. Layers are specifed by
!             histdep(i).
!            
!          time distribution. At each web sector. Layers are the
!              same as the lateral hist.
!               
!       

      subroutine xBgRun
      implicit none
#include "Zmaxdef.h"
#include "Zglobalc.h"
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
#include "Zprivate.h"
#include "Zprivate1.h"
#include "Zprivate2.h"
#include "Zprivate3.h"     
#include "Zprivate4.h" 
#if defined (DOMPI)
#include "mpif.h"
#include "Zmpi.h"
      real*8  ager(maxmpisize)
      real*8  esizer(maxmpisize) 
      real*8  wage, wsize
      real*8  etime1, etime2 
      character*10 numbr
      integer lengnr, rank, rankc
!         to gatter sum of the  all rank data 
      real nrfaiRecA(nrbin, nfai, 4, nsites)
      real nrfaiAllA(nrbin, nfai, 4, nsites)
      real dErfaiA(nrbin, nfai,  nsites)
      real*8 NgA(MaxNoOfAsSites)
      real*8 NeA(MaxNoOfAsSites)
      real*8 NmuA(MaxNoOfAsSites)
      real*8 NhadA(MaxNoOfAsSites)
      real*8 SumElossA(MaxNoOfAsSites)
#define  REALLIMITg  15000
#define  REALLIMITe  9000
#define  REALLIMITmu 7500
#define  REALLIMITh  7500
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
      integer i, j, k, icon, ir, l, m
      integer NN
      integer iij, code, codex
      integer i1, i2, ic, ifai, ji

      real*8 r, Eloss, rinmu, cosang
      real*8 webmin(nrbin, nfai, nsites)
      real*8 dedt, rho, dist, disto, dedxmu
      real*8 fai, rminm

      real*8 u
      logical accept

      real*8 wx, wy, wz
      real   za
      real  de, Ek, f, molu
      real*8  cvh2den
      integer ldep, ridx, faiidx, depidx
      integer lengenv
      integer  ncpu, mcpu ! no. of smashed skeletons, and actully used skeletonsn
      integer  margin
      real*4 wgt, wwgt
      real*4 enhance0      ! since we use only mcpu, the result must be enhanced 
                          ! by a factor of ncpu/mcpu
      real*4 enhance      ! enhance0 * wgt
      integer binw
      character*9 ptcl(4)
      data ptcl/"Photons", "Electrons","Muons", "hadrons"/
      character*9 ptcl2(3)
      data ptcl2/"Electrons", "Muons","All"/
      real power(4)
      integer nstr
      real E0,  cosz, limit(4), limite(4)
      data power/1.,1.,1.,1./
      real  power2(4)
      data power2/1.,1.,1.,1./
      character*128 title
      character*96 evid(nsites), plotid
      real*8 cog, cog2, sumne,  obstimes, Savederg(5)
      real*8 firstz, dd
      real*8 Fai0
      real*8 getFai
      logical dosort
      real prob, tmin, dt
!     ***********************
#include "interface2.h"
!     *********************
!     example
!       histdep:  2 5 6 7 10 /
!       depth   1000 2000 3000 4000 5000 6000 7000 8000 9000 10000
!       ansites = 5
!       w2hl: 0 1 0 0 2 3 4 0 0 5  0 0 ...
      


      lengenv = kgetenv2("LIMIT", input) 
      read(input(1:lengenv),*)  limit
      lengenv = kgetenv2("BINW", input) 
      read(input(1:lengenv),*)  binw
      lengenv = kgetenv2("NCPU", input) 
      read(input(1:lengenv),*)  ncpu
      lengenv = kgetenv2("MCPU", input)
      read(input(1:lengenv),*)  mcpu
      lengenv = kgetenv2("MARGIN", input)
      read(input(1:lengenv),*)  margin
      lengenv = kgetenv2("SeeLowdE", input)
      if(lengenv .le. 0)  then
         write(0,*) ' SeeLowdE not given'
         stop
      endif
      SeeLowdE = input(1:lengenv) .eq. "yes"
      lengenv = kgetenv2("KeepWeight", input)
      KeepWeight = input(1:lengenv) .eq. "yes"

      enhance0 = ncpu
      enhance0 = enhance0/mcpu
      do i = 1, 4
         limite(i) = limit(i) * enhance0 
               !     this must be enhanced from the firs.
      enddo


      write(0,*) ' *** ncpu=',ncpu,
     *            ' mcpu=',mcpu, ' enhance factor=',enhance0
      write(0,*) ' margin=',margin
      do i = 1, NoOfSites
         w2hl(i) = 0
         w2il(i) = 0
      enddo
      do i = 1, NoOfSites
         if(histdep(i) .eq. 0) exit
         hnsites = i
         w2hl( histdep(i)) = i
      enddo

      do i = 1, NoOfSites
         if( indivdep(i) .eq. 0) exit
         ansites = i
         w2il( indivdep(i) ) = i
      enddo

!        now hnsites is the  actual number of sites where statistics
!        by histogram is needed.
!        mapping; depth index to array index; where2loc
!     

 
      call kwhistso( binw )  ! specify binary or ascii write of histogaram
                             !  1--> ascii  2--> binary


      rminm = rmin/10.**(bin/2.0d0)
      r=rmin
      dr = 10.**bin 
      rbin(1) = 0.
      do i = 2, nrbin
         rbin(i) = r
         r = r* dr
      enddo

      return      
!    ******************
      entry xihist
!
!

!     histogram: instanciate
!         rspec (lateral):  
!       

      if(tkarspec) then
         do i = 1, hnsites
            do k = 1, 4
               do j = 1, nfai
                  call kwhisti(rspec(k, j,  i),
     *                 rmin, bin, nrbin,  b'00011' )
                  write(plotid, 
     *            '("Lateral dist. of ",a," at given fai")')
     *                 ptcl(k)
                  call kwhistai(rspec(k, j, i), 
     *            plotid,  "ar", "ptcls", .true.,  power(k), 
     *            "r", "m.u")
               enddo
            enddo
         enddo
      endif
         
      if(tkrtspec) then
          do i = 1, hnsites
             l = histdep(i)
             call cminTime2WebSec(ObsSites(l).pos.xyz,
     *       l, i,   webmin)
          enddo

          do i = 1, hnsites
             do j = 1, 4
!                at  center
                call kwhisti( tspec0(j, i),
     *               -5., 0.1, 200, b'00000')
                call kwhistai( tspec0(j, i),
     *          "Arrival time dist. of "//ptcl(j)//" at center",
     *          "t", "ptcls", .false., 0., "time", "ns")

                do ir=1, nrbin
                   do ifai=1, nfai
                      tmin = webmin(ir, ifai,i)
                      dt = 0.01*10.0**(bin*(ir-1))*100.  ! approx core distnace m
                      dt = dt**0.675*1.e9/3.0e8/30.     ! if sqrt 1m-->0.1 ns 10 m-->0.3 ns
                                       ! 100m 1ns 1km 3ns   4km 6ns
                                       ! dt**0.675  makes larger bin at large distance (<=x2)
                      if(j .eq. 4) dt=dt*10.0 ! for delayed hadrons
                      dt= max(dt, 0.1)
                      call kwhisti( tspec(j, ir, ifai, i),
     *                tmin,  dt, 1000,   b'00000')
                      call kwhistai( tspec(j, ir, ifai,  i),
     *              "Arrival time of "//ptcl(j)//" at (r,faI)",
     *              "rt", "ptcls", .false., 0., "time", "ns")
                  enddo
               enddo
            enddo
         enddo
      endif
      if(tkrtspec) then
         do i = 1, hnsites
            do j = 1, 4
               call kwhistc(tspec0(j,i))
               do ir = 1, nrbin
                  do ifai = 1, nfai
                     call kwhistc(tspec(j, ir, ifai, i))
                  enddo
               enddo
            enddo
         enddo
      endif

      if(tkarspec) then
         do i = 1, hnsites
            do j = 1, 4 
               do ifai = 1, nfai
                  call kwhistc(rspec(j, ifai, i))
               enddo
            enddo
         enddo
      endif
      return
!    ******************
      entry xclearNrfai
!
      
      do i = 1, NoOfSites
         do j= 1, 4
            do k = 1, nfai
               do l = 1,  nrbin
                  nrfaiAll(l, k, j, i) = 0.
                  dErfai(l, k, i) = 0.   !  no j (for all charged ptcls)  > 500 keV
                  if(SeeLowdE) dErfai2(l, k, i) = 0.   !  no j (for all charged ptcls) < 500 keV
               enddo
            enddo
         enddo
      enddo

      do i = 1, ansites
         do j= 1, 4
            do k = 1, nfai
               do l = 1,  nrbin
                  nrfaiRec(l, k, j, i) = 0.
               enddo
            enddo
         enddo
      enddo


      return
!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      entry xBgEvent
#if defined (DOMPI)
      etime1 = MPI_WTIME()
#endif      

      call cqIncident(inci, angle)


#if  FNODATDEF > 0
      bufc=0 
#endif

      cosz = -angle.r(3)
      Fai0 = getFai( -angle.r(1), -angle.r(2)) ! in deg.
!          to rotate coordinate and direction cos in the
!          Obsplane =1 system (x is directed to magnetic east
!          y is to the magnetic  north) to the one in the 
!          system with the x-axis directed to the incident
!          direction. All outuput from this program is
!          measured in this system.  (web data too).
!          Z* = Z exp(-(Fai0)) ;  Z in the detector system
!          Z* in the web system.
!
      CosRot = cos(Fai0*Torad)
      SinRot = sin(Fai0*Torad)
!       
!
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
!          here we use enhanced limit
      call howmuch(limite, E0, NN, cosz)

!      
      do i = 1, NoOfSites
         SumEloss(i) = 0.
         Ng(i) = 0.
         Ne(i) = 0.
         Nmu(i) = 0.
         Nhad(i) = 0. 
      enddo

      obstimes = 0.

      return
!     ***************
      entry xObs(aTrack, id)
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
!    ------------------------------------------
      if(Fai0 .ne. 0. ) then
!         rotate to the web system
         call det2web( aTrack.vec.w.r(1), aTrack.vec.w.r(2),
     *                 aTrack.vec.w.r(1), aTrack.vec.w.r(2))
         call det2web( aTrack.pos.xyz.x, aTrack.pos.xyz.y,
     *       aTrack.pos.xyz.x, aTrack.pos.xyz.y)
      endif
!    ------------------------------------------
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
#if defined (DOMPI)
         write(0, *)'rank=',mpirank,'  obstimes=', obstimes,
     *    ' ptclE=',aTrack.p.fm.p(4)
#else
         write(0, *),'  obstimes=', obstimes,
     *    ' ptclE=',aTrack.p.fm.p(4)
#endif
         do i = 1, min(4,Stack_pos)
            write(0,*)' stack tops=', Stack(i).p.fm.p(4)
         enddo
      endif
!     ***************
      code = aTrack.p.code
      ldep =  aTrack.where

      if(id .eq. 2 .and. code .le. 6) then
         codex=min(code, 4)
         wz = aTrack.vec.w.r(3)  ! downgoing < 0
         if(wz .gt. 0) return
         wz = -wz
         if( KeepWeight ) then
            wgt = aTrack.wgt
            enhance= enhance0*wgt
         else
            wgt = 1.
            enhance= enhance0
         endif

         if(code .eq. kphoton) then
            Ng(ldep) = Ng(ldep) + enhance
         elseif(code .eq. kelec) then
            Ne(ldep) = Ne(ldep) + enhance
         elseif(code .eq. kmuon) then
            Nmu(ldep) = Nmu(ldep) + enhance
         else
            Nhad(ldep)=Nhad(ldep) + enhance
         endif

         r = sqrt( aTrack.pos.xyz.x**2 +
     *                 aTrack.pos.xyz.y**2 )

         Ek = aTrack.p.fm.p(4) -aTrack.p.mass

!       ------------- compute energy loss rate
         if(aTrack.p.charge .ne. 0  ) then

!             -----------------
!                      /|    |
!                     / |   1g/cm2k
!                    /A |    |
!            -------------------
!                  / ptcl direction  
!         get energy loss when aTrack goes some distance
!         of which vertical gramage is 1g/cm2.
!         Gramage the particle travel is 
!         1/cos where cos is the cos of angle (i.e, A if Fig)
!          in the detctor system.
!         1g/cm^2 = 10-3kg/10-4 m^2 =10 kg/m^2.
!         To travel  1 g/cm^2  along shower axis, the ptcl must
!         run dist kg/m^2
            if(abs(wz) .gt. 1.d-2) then
               dist =10./wz    ! in kg/m2/(g/cm2)
            else
!                 for safety
               dist =1000.
            endif
            if( HEobs ) then
!                   the ptcls is the one obsrved at skeleton making time
!                   we must compute dE/dt here
               rho = cvh2den(aTrack.pos.height)
               call cdedxInAir(aTrack.p, rho, dedt) ! dedt; GeV/(kg/m2)
               if(aTrack.p.code .eq. kmuon ) then
!                 dE/dx due to muon pair, brem, nuc.i
                  call cmudEdx(MuNI, MuBr, MuPr, aTrack.p.fm.p(4),
     *                 dedxmu)  ! dedxmu in GeV/(g/cm2)
                  dedxmu = dedxmu /10. !  GeV/(kg/m2)
                  dedt = dedt + dedxmu
               endif
            else
!               we can use already computed one
               call cqElossRate(dedt) !  loss rate GeV/(kg/m^2)
            endif
!                       energy in 1 g/cm2 of vertical direction
            Eloss =min( dedt*dist, dble(Ek)) * enhance  !  GeV/(g/cm2)

            SumEloss(ldep)=SumEloss(ldep) + Eloss
         else
            Eloss=0.
         endif
!          
!          get web sector indexes (r and fai index )
!
         molu = ObsSites(ldep).mu
         rinmu =r/molu
         sr = rinmu
         if(rinmu .gt. rminm) then
            ridx= log10(rinmu/rminm)/bin + 1
         else
            ridx =0
         endif
!              fai is  in    -15 to 345  (for dfai=30.)
         fai=getFai(aTrack.pos.xyz.x, aTrack.pos.xyz.y)
         fai= mod(fai + 360.d0,   360.d0)
         if(fai .gt. (360.d0-dfai/2.0d0)) fai= fai-360.d0
         faiidx=(fai+dfai/2.0d0) /dfai + 1
!
         if(ridx .gt. 0 .and. ridx .le. nrbin ) then
!               for all particles
            nrfaiAll(ridx, faiidx, codex, ldep) =
     *           nrfaiAll(ridx, faiidx, codex, ldep) + enhance
            if(SeeLowdE) then 
!                 separate low energy contribution ; Eloss = 0 for neutral 
               if(Ek .gt. 500.d-6) then
                  dErfai(ridx, faiidx, ldep) =
     *                 dErfai(ridx, faiidx,  ldep) +  Eloss
               else
                  dErfai2(ridx, faiidx, ldep) =
     *                 dErfai2(ridx, faiidx,  ldep) +  Eloss
               endif
            else
!                don't separate low E and high E ptcls
               dErfai(ridx, faiidx, ldep) =
     *              dErfai(ridx, faiidx,  ldep) +  Eloss
            endif  
!
!             ==============
!               for  individual particle
            ji = w2il(ldep)     !  indivdep index.
            if(ji .gt.  0 ) then
               prob =  recprob(ridx, codex, ji)
               if(prob .gt. 1.) then
                  wwgt = wgt
                  accept = .true.
               else
!                 Example:  if prob=10^-3 and wgt=2x10^3, 
!                  record it with weight =2 (=prob*wgt)
                  prob = prob*wgt ! wgt =1 or wgt > 1
                  if(prob .gt. 1.) then
                     wwgt = prob
                     accept = .true.
                  else
                     call rndc(u)
                     if(u .lt. prob) then
                        accept =.true.
                        wwgt=1.
                     else
                        accept = .false.
                     endif
                  endif                  
               endif  ! computtion of accept end
               if(accept) then
                  nrfaiRec(ridx, faiidx, codex, ji) =
     *              nrfaiRec(ridx, faiidx, codex, ji) + wwgt
#if  FNODATDEF > 0
                 if(bufc .lt. bufsize) then
                    bufc = bufc + 1
                    buf(bufc).ldep=ldep
                    buf(bufc).code=code
                    buf(bufc).subcode=aTrack.p.subcode
                    buf(bufc).charge = aTrack.p.charge
                    buf(bufc).ridx=ridx
                    buf(bufc).faiidx= faiidx
                    buf(bufc).rinmu = rinmu
                    buf(bufc).fai= fai
                    buf(bufc).Ek = Ek
                    buf(bufc).t = aTrack.t
                    buf(bufc).wx=-aTrack.vec.w.r(1)
                    buf(bufc).wy=-aTrack.vec.w.r(2)
                    buf(bufc).wz=wz
                    buf(bufc).wgt = wwgt !  not wgt
                 else
#if BUFSIZE > 200000
                    write(0,*) ' too many partciles --> buf'
                    write(0,*) ' you must make BUFSIZE in interface.f'
                    write(0,*) ' to the standard value (<2x10^5) '
                    stop
#else
                    write(fnodat) bufc, buf
                    bufc= 0
#endif
                 endif  
#else
                 write(*,
     *  '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p,2f8.4,f10.6,1pE11.3)'
     *    )
     *          ldep,  code,  aTrack.p.subcode,
     *          aTrack.p.charge, ridx, faiidx,
     *          rinmu, fai,
     *          Ek, aTrack.t, 
     *          -aTrack.vec.w.r(1), 
     *          -aTrack.vec.w.r(2),  wz, wwgt  ! not wgt
#endif
               endif  ! end accepted ptcl treatment
            endif   ! end individual ptcl output treatment
         endif  ! ridx  >0 ended



!        ================ for histograming
         i=w2hl(ldep)
         if(i .gt. 0 ) then
            if(tkarspec) then
               call kwhist( rspec(codex,faiidx, i),  sr, enhance)
            endif

            if( tkrtspec ) then
               if( ridx .eq. 0 ) then
                  call kwhist( tspec0(codex, i),
     *              sngl( aTrack.t ), enhance)
               elseif(ridx .le. nrbin) then
                  call kwhist( tspec(codex, ridx, faiidx,  i),
     *                 sngl( aTrack.t ), enhance)
               endif
            endif
         endif       
      endif   ! code and id judge
      return
!     **************
      entry xEnEvent
!     **************
#if FNODATDEF > 0
#if BUFSIZE < 200000 
      if(bufc .gt. 0) then
         write(fnodat)  bufc, buf
         bufc=0
      endif
      close(fnodat)
#endif
#endif


      firstz= Zfirst.pos.depth*0.1

#if defined (DOMPI)
      write(0,*) ' rank=',mpirank, ' closed main file'
      etime2 = MPI_WTIME()
      write(0,*) ' elapsed time for this event = ',
     * etime2-etime1
      call MPI_Barrier(MPI_COMM_WORLD, icon)
#include "inc_gatherNrfai.f"
#include "inc_gatherHyb.f"
      if( mpirank .eq. 0) then
         write(0,*) " rank 0 is reading all data and combining now"
#include "inc_readAndput.f"
      endif  
#endif
#if defined (DOMPI)
      if( mpirank .eq.  0) then
#endif
#include "inc_writeHyb.f"
#include "inc_writeNrfai.f"
#if defined (DOMPI)
      endif
#endif

!       fwollowings are not supported in MPI job.
#if ! defined (DOMPI)
      do i = 1, hnsites
         j = histdep(i)
         write( evid(i), 
     *   '(i3, f7.1,  f6.3, f6.3, i5,  i4)')
     *   histdep(i), ASDepthList(j)*0.1,
     *   ASObsSites(j).age, ASDepthList(j)*0.1/cog,
     *   int(ASObsSites(j).mu), int(cog)
      enddo

      if(tkarspec) then
         do i = 1, hnsites
            k=histdep(i)
            do j = 1, 4
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( DepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               do l = 1, nfai
                  call kwhistdir(rspec(j, l,  i),  dirstr)
                  call kwhists(  rspec(j, l, i), 0. )
                  call kwhistev( rspec(j, l, i), EventNo)
                  call kwhistid( rspec(j, l,  i), evid(i))
                  call kwhistp( rspec(j, l, i),  fno)
!                        *********** deallocate ********           
                  call kwhistd( rspec(j, l, i) )
               enddo  ! code loop
            enddo ! fai loop
         enddo ! depth loop
      endif

      if( tkrtspec ) then
         do i = 1, hnsites
            do j = 1, 4   
               call kwhists( tspec0(j,i), 0.)
               call kwhistev( tspec0(j,i), EventNo)
               call kwhistid( tspec0(j,i), evid(i))
               k=histdep(i)
               dirstr = " "
               write( dirstr,'(a,"/d",i4, "/")')
     *              ptcl(j), int( ASDepthList(k)*0.1 )
               call kseblk( dirstr, "|", nstr)
               call kwhistdir( tspec0(j,i),  dirstr )
               call kwhistp( tspec0(j,i),  fno )
!                 *********** deallocate ********                  
               call kwhistd( tspec0(j,i) )

               do ifai= 1, nfai
                  do ir= 1, nrbin
                     call kwhists( tspec(j,ir, ifai,i), 0.)
                     call kwhistev(tspec(j,ir, ifai,i), EventNo)
                     call kwhistid( tspec(j,ir, ifai,i), evid(i))
                     dirstr = " "
                     write(dirstr,'(a,"/d",i4, "/F",i2,"/")')
     *                    ptcl(j), int( DepthList(k)*0.1), ifai
                     call kseblk(dirstr, "|", nstr)
                     call kwhistdir(tspec(j,ir, ifai,i),  dirstr)
                     call kwhistp( tspec(j,ir, ifai,i),  fno)
!                        *********** deallocate ********           
                     call kwhistd( tspec(j, ir, ifai, i) )
                  enddo
               enddo
            enddo  ! code loop
         enddo   ! depth loop
      endif
#endif

      end
      subroutine watchdog(h)
      implicit none
#include "Zmaxdef.h"
#include "Zobs.h"
#include "Zprivate.h"
#include "Zprivate2.h"
      type(histogram1):: h
      integer i
      write(0,*) ' h%c%init=', h%c%init
      write(0,*) ' h%c%title=', h%c%title
      if( h%c%init .eq. 'initend' ) then
         do i = 1, 20
            write(0,*) 'i=',i, ' h%dnw(i)=',h%dnw(i)
         enddo
      endif
      end
      real*8 function getFai(x, y)
      implicit none
#include "Zglobalc.h"
!        get polar angle of the point (x,y) in deg.
      real*8 x, y  ! input postion or direction-cos components
      getFai = atan2(y, x)*Todeg
      end
      subroutine det2web(x, y, xs, ys)
      implicit none 
#include "Zmaxdef.h"
#include "Zobs.h"
#include "Zprivate.h"
      real*8 x, y ! input.  x,y component of the coordinate
                  !   or direction cos  in the detctor system.
      real*8 xs,ys ! ouput. values transformed to web system
               !  where the x-axis is directed to the incident
               !  xs,ys can be x,y
      
!      Z* = Z * exp(-Fai0) = (x + iy)(cos-i sin)
!                          = x cos  + y sin + i ( y cos -x sin)
!
      real*8 temp
      temp = x*CosRot + y*SinRot
      ys = y*CosRot - x*SinRot
      xs = temp
      end
