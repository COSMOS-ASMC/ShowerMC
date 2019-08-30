!  FNODATDEF   if undef,  ascii main output; else  binary output:
!        execSSHtemplate.sh   
!   or   execSGEtemplate.sh must be modified.
!
#define FNODATDEF 33
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
#include "Zprivate.h"
#include "Zprivate1.h"
#include "Zprivate2.h"
#include "Zprivate3.h"     

      type(track):: incident
      type(coord):: AngleAtObs

      logical  HEobs            ! if T, currently observing 
      common /ZHEobs/ HEobs     !  particle is the one  obsrved at skeelton making time

            

            
      save rspec, lossrspec, arspec,  respec
      save rzspec,  zfspec, rtspec1, rtspec2, retspec1, retspec2
      save rezspec,  rzfspec, rfspec, efspec, refspec


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
      integer  ncpu, mcpu ! no. of smashed skeletons, and actully used skeletonsn
      integer  margin
      real*4 enhance      ! since we use only mcpu, the result must be enhanced 
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
      save


!     
!     ***********************
#include "interface2.f"
!     *********************
!     example
!       histdep:  2 5 6 7 10 /
!       depth   1000 2000 3000 4000 5000 6000 7000 8000 9000 10000
!       ansites = 5
!       w2hl: 0 1 0 0 2 3 4 0 0 5  0 0 ...
      

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
      do i = 1, nsites
         w2hl(i) = 0
         w2il(i) = 0
      enddo
      do i = 1, nsites
         if(histdep(i) .eq. 0) exit
         hnsites = i
         w2hl( histdep(i)) = i
      enddo

      do i = 1, nsites
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
      entry ihist
!
!

!     histogram: instanciate
!         rspec (lateral):  
      if(tklat) then
         do i = 1, hnsites
            do j = 1, 4         ! g,e,mu,had
               call kwhisti(rspec(j, i),
     *              rmin, bin, nrbin,  b'00011' )
               call kwhistai(rspec(j,i), 
     *         "Lateral Dist. of "//ptcl(j),
     *         "lat", "ptcls", .true., power(j), "r", "m.u")
            enddo
         enddo
      endif


!        elosrspec (energy loss lateral)  10m-10km log bin 0.1. e+mu, e,mu
!          
!       
      if(tkelosslat) then
         do i = 1, hnsites
            do j = 1, 3         !   e, mu, e+mu
               call kwhisti( lossrspec(j, i),
     *          rmin, bin, nrbin, b'00011')
               call kwhistai(lossrspec(j, i), 
     *         "dE/dx lateral dist. of "//ptcl2(j),
     *         "dEdxLat", "GeV/(g/cm^2)", .true., power2(j),
     *         "r", "m.u")
            enddo
         enddo   
      endif


      if(tkarspec) then
         do i = 1, hnsites
            do j=1, 4
               call kwhisti2(arspec(j, i),
     *              0.,  30.0, 12,     b'00010',
     *              rmin, bin, nrbin,  b'00011' )
               call kwhistai2(arspec(j,i), 
     *         "Lateral dits. of "//ptcl(j)//" with given Fai bin",
     *         "ar", "ptcls", .true.,  power(j), 
     *         "azimuth", "deg", "r", "m.u")
            enddo
         enddo
      endif
         

!
      if(tkrespec) then
         do i = 1, hnsites
            do j= 1, 2
               call kwhisti2( respec(j, i),
     *             0.01, 0.2, 20,    b'00011', 
     *             500.e-6, 0.1,    50,    b'00001')
               call kwhiststep2(respec(j, i), 2)
               call kwhistai2( respec(j, i),
     *         "Energy Spec. of "//ptcl(j)//" at  diff. r",
     *         "re", "ptcls", .true.,  1., 
     *         "r", "m.u", "E", "GeV")
            enddo
!                 mu            
            call kwhisti2( respec(3, i),
     *            0.01, 0.2, 20       , b'00011', 
     *            0.031627,  0.1,  38, b'00011' )
            call kwhiststep2(respec(3, i), 2)
            call kwhistai2( respec(3, i),
     *      "Energy Spec. of mu at diff. r",
     *      "re", "ptcls", .true.,  0.,
     *      "r",  "m.u", "E", "GeV")
         enddo
      endif


      if(tkrzspec) then
         do i = 1, hnsites
            do j= 1, 2
               call kwhisti2( rzspec(j, i),
     *           0.01,   0.2,  20,  b'00011', 
     *           0.,     1.0,     20,  b'10000')
               call kwhiststep2(rzspec(j, i), 2)
               call kwhistai2( rzspec(j, i),
     *         "Zenith angle dist. of "//ptcl(j)//" at diff. r",
     *         "rz", "ptcls", .true., 0.,
     *         "r", "m.u", "cosz", " ")
            enddo
            call kwhisti2(rzspec(3, i),
     *           0.01, 0.2, 20, b'00011', 
     *           0.,  1.0,   20, b'10000' )  
           call kwhiststep2(rzspec(3, i), 2)
            call kwhistai2( rzspec(3, i),
     *         "Zenith angle dist. of m  with diff. r",
     *         "rz", "ptcls", .true., 0.,
     *         "r", "m.u", "cosz", " ")
         enddo
      endif         


      if(tkzfspec) then
         do i = 1, hnsites
            do j= 1, 3
               call kwhisti2( zfspec(j, i),
     *              0.0,  1.0, 10,  b'10000', 
     *              -1.0, 1.0, 50,  b'10000')
               call kwhiststep2(zfspec(j, i), 2)
               call kwhistai2( zfspec(j, i),
     *         "f=(wx,wy)*(x,y) spectrum of "//ptcl(j)//
     *         " with diff. cosz", 
     *         "zf", "ptcls", .true., 0.,
     *         "cosz", " ", "f", " ")
            enddo
         enddo
      endif         


      if(tkrfspec) then
         do i = 1, hnsites
            do j= 1, 3
               call kwhisti2( rfspec(j, i),
     *                0.01,  0.2, 20,  b'00011', 
     *               -1.0,  1.0, 50,   b'10000')
               call kwhiststep2(rfspec(j, i), 2)
               call kwhistai2( rfspec(j, i),
     *         "f spectrum of "//ptcl(j)//" with diff. r",
     *         "rf", "ptcls", .true., 0.,
     *          "r", "m.u", "f", " ")
            enddo
         enddo
      endif         


      if(tkefspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhisti2( efspec(j, i),
     *              500.e-6,  0.2, 20,  b'00001', 
     *              -1.0,     1.0, 50,  b'10000')
               call kwhiststep2(efspec(j, i), 2)
               call kwhistai2( efspec(j, i),
     *         "f spectrum of "//ptcl(j)//" with diff. E",
     *         "ef", "ptcls",.true., 0.,
     *         "E", "GeV", "f", " ")
            enddo
         enddo
      endif         

      if(tkrtspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhisti2( rtspec1(j, i),
     *              1.0, 0.2,  10,   b'00011', 
     *              10.,  0.1, 55, b'00011'  )
               call kwhiststep2(rtspec1(j, i), 2)
               call kwhistai2( rtspec1(j, i),
     *         "Arrival time dist. of "//ptcl(j)//" at diff. r",
     *         "rt", "ptcls", .true., 0.,
     *         "r", "m.u", "time", "ns")

               call kwhisti2( rtspec2(j, i),
     *             0.,  1.0, 10,  b'10000', 
     *             0., 0.25, 300, b'00000' )
               call kwhiststep2(rtspec2(j, i), 2)
               call kwhistai2( rtspec2(j, i),
     *         "Arrival time dist. of "//ptcl(j)//" at diff. r",
     *         "rt", "ptcls", .true., 0.,
     *         "r", "m.u", "time", "ns")
            enddo
         enddo
      endif

      if(tkretspec)  then
         do i = 1, hnsites
            do j = 1, 2
!                  g,e
               call kwhisti3(retspec1(j, i), 
     *          1.0,     0.2,  10,    b'00011',
     *          500.e-6, 0.25, 15,    b'01001',
     *          10.0, 0.1,  55,      b'00011' )
               call kwhiststep3(retspec1(j, i), 2, 2)
               call kwhistai3(retspec1(j, i), 
     *         "Arrival time dist. of "//ptcl(j)//" with diff. r&E",
     *         "ret", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "Time", "ns")

               call kwhisti3(retspec2(j, i), 
     *          0.0,  1.0,     10,   b'10000',
     *          500.e-6, 0.25, 12,  b'01001',
     *          0.0,  0.25,   200,  b'00000' )
               call kwhiststep3(retspec2(j, i), 2,2)
               call kwhistai3(retspec2(j, i), 
     *         "Arrival time dist. of "//ptcl(j)//" with diff. r&E",
     *         "ret2", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "Time", "ns")
            enddo
!               mu
            call kwhisti3(retspec1(3, i), 
     *           1.0,   0.2,    10,       b'00011',
     *           100.e-3, 0.25, 10,       b'01011',
     *           10.0,    0.1,  72,       b'00011' )
            call kwhiststep3(retspec1(3, i),  2, 2)
            call kwhistai3(retspec1(3, i), 
     *         "Arrival time dist. of "//ptcl(3)//" with diff. r&E",
     *         "ret", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "Time", "ns")


            call kwhisti3(retspec2(3, i), 
     *           0.0,  1.0,  10,    b'10000',
     *           100.e-3,    0.25, 10,   b'01011',
     *           0.0,  0.25,    200,     b'00000' )
            call kwhiststep3(retspec2(3, i), 2,2)
            call kwhistai3(retspec2(3, i), 
     *         "Arrival time dist. of "//ptcl(3)//" with diff. r&E",
     *         "ret2", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "Time", "ns")
         enddo
      endif


      if(tkrezspec) then
         do i = 1, hnsites
            do j = 1, 2
               call kwhisti3( rezspec(j, i), 
     *         0.1,   0.2,  15,     b'00011',
     *         500.e-6, 0.25, 14,   b'01001', 
     *         0.0, 1.0,   20,    b'10000')
               call kwhiststep3(rezspec(j, i), 3,2)
               call kwhistai3(rezspec(j, i), 
     *         "cos zenith dist. of "//ptcl(j)//" with diff. r&E",
     *         "rez", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "cosz", " ") 
            enddo
            call kwhisti3(rezspec(3, i), 
     *         0.1,  0.2,   15,   b'00011',
     *         100.e-3, 0.25, 10, b'01011', 
     *         0.0, 1.0,     20,   b'10000')
            call kwhiststep3(rezspec(3, i), 3,2)
            call kwhistai3(rezspec(3, i), 
     *         "cos zenith dist. of "//ptcl(3)//" with diff. r&E",
     *         "rez", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "cosz", " ") 
         enddo
      endif



      if(tkrzfspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhisti3(rzfspec(j, i), 
     *         0.1, 0.2,        15,   b'00011',
     *         0.0, 1.0,        10,   b'10000',  
     *         -1.0, 1.0,       20,   b'10000')
               call kwhiststep3(rzfspec(j, i), 3,2)
               call kwhistai3(rzfspec(j, i), 
     *         "f spectrum of "//ptcl(j)//" with diff r&cosz",
     *         "rzf", "ptcls", .true., 0.,
     *         "r", "m.u", "cosz", " ", "f", " ")
            enddo
         enddo
      endif

      if(tkrefspec) then
         do i = 1, hnsites
            do j = 1, 2
               call kwhisti3( refspec(j, i), 
     *         0.1,  0.2,      15,   b'00011',
     *         500.e-6, 0.25,  16,    b'01001',  
     *         -1.0, 1.0,     20,   b'10000')
               call kwhiststep3(refspec(j, i), 3,3)
               call kwhistai3(refspec(j, i),
     *         "f spectrum of "//ptcl(j)//" with diff. r&E",
     *         "ref", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "f", " ")
            enddo
            call kwhisti3(refspec(3, i), 
     *       0.1, 0.2,      15,    b'00011',
     *       100.e-3,  0.25, 10,    b'01001',  
     *        -1.0, 1.0,     20,   b'10000')
            call kwhiststep3(refspec(3, i), 3,2)
            call kwhistai3(refspec(3, i),
     *         "f spectrum of "//ptcl(3)//" with diff. r&E",
     *         "ref", "ptcls", .true., 0.,
     *         "r", "m.u", "E", "GeV", "f", " ")
         enddo
      endif

      return
!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      entry xBgEvent


      call cqIncident(inci, angle)


#if  FNODATDEF > 0
      bufc=0 
#endif

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
!          here we use enhanced limit
      call howmuch(limite, E0, NN, cosz)

!      
      do i = 1, NoOfASSites
         SumEloss(i) = 0.
         Ng(i) = 0.
         Ne(i) = 0.
         Nmu(i) = 0.
         Nhad(i) = 0.
      enddo

      if(tklat) then
         do i = 1, hnsites
            do j = 1, 4
               call kwhistc(rspec(j,i))
            enddo
         enddo
      endif

      if(tkelosslat) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc(lossrspec(j,i))
            enddo
         enddo
      endif


      if(tkarspec) then
         do i = 1, hnsites
            do j = 1, 4
               call kwhistc2(arspec(j,i))
            enddo
         enddo
      endif

      if(tkrespec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc2(respec(j,i))
            enddo
         enddo
      endif


      if(tkrzspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc2(rzspec(j,i))
            enddo
         enddo
      endif

      if(tkrfspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc2(rfspec(j,i))
            enddo
         enddo
      endif



      if(tkefspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc2(efspec(j,i))
            enddo
         enddo
      endif




      if(tkrtspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc2(rtspec1(j,i))
               call kwhistc2(rtspec2(j,i))
            enddo
         enddo
      endif

         
      if(tkretspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc3(retspec1(j,i))
               call kwhistc3(retspec2(j,i))
            enddo
         enddo
      endif



      if(tkrezspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc3(rezspec(j,i))
            enddo
         enddo
      endif



      if(tkrzfspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc3(rzfspec(j,i))
            enddo
         enddo
      endif


      if(tkrefspec) then
         do i = 1, hnsites
            do j = 1, 3
               call kwhistc3(refspec(j,i))
            enddo
         enddo
      endif

      do i = 1, hnsites
         do j= 1, 4
            do k = 1, nfai
               do l = 1,  nrbin
                  nrfaiRec(l, k, j, i) = 0.
                  nrfaiAll(l, k, j, i) = 0.
                  dErfai(l, k, i) = 0.   !  no j (for all charged ptcls) 
               enddo
            enddo
         enddo
      enddo

      obstimes = 0.

      return
!     ***************
      entry xObs(aTrack, id)
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!
!     **************************
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
      if(id .eq. 2 .and. code .le. 6) then
         codex=min(code, 4)
         wz = aTrack.vec.w.r(3)  ! downgoing < 0
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

         r = sqrt( aTrack.pos.xyz.x**2 +
     *                 aTrack.pos.xyz.y**2 )

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
               dist =10./wz    ! in kg/m2/(g/cm2)
            else
!                 for safety
               dist =1000.
            endif
            if( HEobs ) then
!                   the ptcls is the one obsrved at skeleton making time
!                   we must compute dedt here
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
!                       energy loss rate
            Eloss = dedt*dist*enhance   !  GeV/(g/cm2)

            SumEloss(ldep)=SumEloss(ldep) + Eloss
         else
            Eloss=0.
         endif

         j = w2il(ldep) 
         if(j .eq. 0) goto 100

         molu =ObsSites(ldep).mu
         rinmu =r/molu
         sr = rinmu
         if(rinmu .gt. rminm) then
            ridx= log10(rinmu/rminm)/bin + 1
         else
            ridx =0
         endif
!           fai is  in    -15 to 345  (for dfai=30.)
         fai=atan2(aTrack.pos.xyz.y, aTrack.pos.xyz.x)*
     *        57.29577951308230d0
         fai= mod(fai + 360.d0,   360.d0)
         if(fai .gt. (360.d0-dfai/2.0d0)) fai= fai-360.d0
         faiidx=(fai+dfai/2.0d0) /dfai + 1

         if(ridx .gt. 0 .and. ridx .le. nrbin ) then
            call rndc(u)
            accept = u .le. recprob(ridx, codex, j)
            if(accept) then
               nrfaiRec(ridx, faiidx, codex, j) =
     *              nrfaiRec(ridx, faiidx, codex, j) + 1.0
            endif
            nrfaiAll(ridx, faiidx, codex, j) =
     *           nrfaiAll(ridx, faiidx, codex, j) + 1.0
            dErfai(ridx, faiidx, j) =
     *           dErfai(ridx, faiidx,  j) +  Eloss
         else
            accept = .false.
         endif
         if(accept) then
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
            else
               write(fnodat) bufc, buf
               bufc= 0
            endif
#else
            if(wz .lt. 0.999) then
               write(*,
     *          '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6)')
     *              ldep,  code,  aTrack.p.subcode,
     *              aTrack.p.charge, ridx, faiidx,
     *              rinmu, fai,
     *              Ek, aTrack.t, 
     *              -aTrack.vec.w.r(1),  -aTrack.vec.w.r(2),  wz
            else
               write(*,'(6i3,1pE11.3, 0p, f6.1, 1p2E11.3, a)')
     *              ldep,  code,  aTrack.p.subcode,
     *              aTrack.p.charge, ridx, faiidx,
     *              rinmu, fai,
     *              Ek, aTrack.t,   " 0 0 1"
            endif
#endif
         endif
      endif
 100  continue
!            ------------------
      i=w2hl(ldep)

      if(id .eq. 2 .and. code .le. 6 .and. i .gt. 0) then
         if( tklat ) then
!                    lateral
            call kwhist(rspec(codex, i), sr, enhance )
         endif
         if( tkelosslat  
     *        .and. aTrack.p.charge .ne. 0 ) then
!c            de = Eloss*enhance
            de = Eloss
            call kwhist(lossrspec(3,i), sr, de)
            if( code .eq.  kelec ) then
!                          by electrons
               call kwhist(lossrspec(1,i), sr, de)
            else
!                      by other charged particles
               call kwhist(lossrspec(2,i), sr, de)
            endif
!                  end of eloss late
         endif


         if(tkarspec) then
!                arspec
            call kwhist2( arspec(codex, i), sngl(fai), sr, enhance)
         endif
      endif
      if(id .eq. 2 .and. code .le. 3 .and. i .gt. 0) then
         if( tkrespec ) then
!                  re spectrum
            call kwhist2( respec(code, i), sr, Ek, enhance)
         endif

         if( tkrzspec ) then
            call kwhist2( rzspec(code, i), sr, 
     *           sngl(wz), enhance )
         endif

         if( wz .lt. 1.0d0 ) then
            temp = 1.d0 - wz**2
            if(temp .gt. 0.) then
               temp = sqrt(temp)
               f = (
     *              aTrack.pos.xyz.x* aTrack.vec.w.r(1) +
     *              aTrack.pos.xyz.y* aTrack.vec.w.r(2) ) /r
     *              /temp             
            else
               f = 1.0
            endif
         else
            f = 1.0
         endif
         if( tkzfspec .and. f .lt. 1.0 ) then
            call kwhist2( zfspec(code, i), 
     *           sngl(wz), f, enhance )
         endif


         if( tkefspec .and. f .lt. 1.0 ) then
            call kwhist2( efspec(code, i), 
     *           Ek, f, enhance )
         endif

         if( tkrfspec .and. f .lt. 1.0 ) then
            call kwhist2( rfspec(code, i), 
     *           sr, f, enhance )
         endif

         if( tkrtspec ) then
            call kwhist2( rtspec1(code, i), sr, 
     *           sngl( aTrack.t ), enhance )
            call kwhist2( rtspec2(code, i), sr, 
     *           sngl( aTrack.t ), enhance )
         endif

         if(tkretspec ) then
            call kwhist3( retspec1(code, i), sr,  Ek, 
     *           sngl( aTrack.t ), enhance)
            call kwhist3( retspec2(code, i), sr,  Ek, 
     *           sngl( aTrack.t ), enhance)
         endif

         if( tkrezspec ) then
            call kwhist3( rezspec(code,i), sr, Ek, 
     *           sngl(wz),  enhance ) 
         endif
         
         if( tkrzfspec ) then
            call kwhist3( rzfspec(code, i),  sr, 
     *           sngl(wz), f, enhance )
         endif

         if( tkrefspec ) then
            call kwhist3( refspec(code, i),  sr, 
     *           Ek, f, enhance)
         endif
      endif       
      return
!     **************
      entry xEnEvent
!     **************
#if FNODATDEF > 0
!      if(fnodat .gt. 0) then
         if(bufc .gt. 0) then
            write(fnodat)  bufc, buf
            bufc=0
         endif
!      endif
#endif


!      call cqFirstID(firstz)
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
         else
            write(*,
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
            else
               write(*, '("t ", i3, 2f7.1,  2f6.3,
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
         else
            write(*,*)
         endif
      endif


!        here we write limit; at reduce process, non enhanced value
!        is needed
      write(fnonrfai,'(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,3i4)')
     *   EventNo, E0, NN, cosz, limit, nrbin, nfai, ansites
      do i = 1, ansites
         do j = 1, 4
            do k = 1, nfai
               l = indivdep(i)
               write(fnonrfai, '("rec",f7.1, 4i4)' )
     *          ASDepthList(l)*0.1, l, i, j, k
               write(fnonrfai, '(1p10E11.3)')
     *             ( nrfaiRec(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo

      do i = 1, ansites
         do j = 1, 4
            do k = 1, nfai
               l = indivdep(i)
               write(fnonrfai, '("all",f7.1, 4i4)' )
     *         ASDepthList(l)*0.1, l, i, j, k
               write(fnonrfai, '(1p10E11.3)')
     *             ( nrfaiAll(l,k,j,i)*enhance, l=1,nrbin )
            enddo
         enddo
      enddo
!           dE/dx lateral
      do i = 1, ansites
         do k = 1, nfai
            l = indivdep(i)
            write(fnonrfai, '("dE/dx",f7.1, 3i4)' )
     *         ASDepthList(l)*0.1, l, i, k
            write(fnonrfai, '(1p10E11.3)')
     *             ( dErfai(l,k,i)*enhance, l=1,nrbin )
         enddo
      enddo

      write(fnonrfai, *)


      do i = 1, hnsites
         j=histdep(i)
         write(evid(i), 
     *   '(i3, f7.1,  f6.3, f6.3,
     *   i5,  i4)')
     *   histdep(i), ASDepthList(j)*0.1,
     *   ASObsSites(j).age, ASDepthList(j)*0.1/cog,
     *   int(ASObsSites(j).mu), int(cog)
      enddo


      if( tklat ) then
         do j = 1, 4
            do i = 1, hnsites
               k = histdep(i)
               call kwhists( rspec(j,i), 0. ) ! 0. means area norm.
               call kwhistev(rspec(j,i), EventNo)
               call kwhistid(rspec(j,i), evid(i))
               call kwhistdir(rspec(j,i), ptcl(j)//"/")
               call kwhistp( rspec(j,i),  fno)
            enddo
         enddo
      endif

      if(tkelosslat) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists( lossrspec(j,i), 0. ) ! area norm
               call kwhistev(lossrspec(j,i), EventNo)
               call kwhistid( lossrspec(j, i), evid(i))
               call kwhistdir( lossrspec(j, i), ptcl2(j)//"/")
               call kwhistp( lossrspec(j, i), fno)
            enddo
         enddo

      endif
 120  continue
      
      if(tkarspec) then
         do i = 1, hnsites
            k=histdep(i)
            do j = 1, 4
               call kwhists2( arspec(j, i), 0. )
               call kwhistev2( arspec(j,i), EventNo)
               call kwhistid2( arspec(j, i), evid(i))
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(arspec(j, i),  dirstr)
               call kwhistp2( arspec(j, i),  fno)
            enddo
         enddo
      endif
      
      if(tkrespec) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( respec(j, i), 0. )
               call kwhistev2( respec(j,i), EventNo)
               call kwhistid2( respec(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(respec(j, i),  dirstr)
               call kwhistp2( respec(j, i),  fno)
            enddo
         enddo
      endif


      if( tkrzspec ) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( rzspec(j, i), 0.)
               call kwhistid2( rzspec(j, i),  evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2( rzspec(j, i),  dirstr)
               call kwhistp2( rzspec(j, i),  fno)
            enddo
         enddo
      endif

      
      if( tkzfspec ) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( zfspec(j, i), 0.)
               call kwhistid2( zfspec(j, i), evid(i) )
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(  zfspec(j, i),  dirstr)
               call kwhistp2( zfspec(j, i),  fno)
            enddo
         enddo
      endif


      
      if( tkrfspec ) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( rfspec(j, i), 0.)
               call kwhistid2( rfspec(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(rfspec(j, i),  dirstr)
               call kwhistp2( rfspec(j, i),  fno)
            enddo
         enddo
      endif

      
      if( tkefspec ) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( efspec(j, i), 0.)
               call kwhistev2(efspec(j,i), EventNo)
               call kwhistid2( efspec(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(efspec(j, i),  dirstr)
               call kwhistp2( efspec(j, i),  fno)
            enddo
         enddo
      endif


      if( tkrtspec ) then
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( rtspec1(j,i), 0.)
               call kwhistev2( rtspec1(j,i), EventNo)
               call kwhistid2( rtspec1(j,i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(rtspec1(j,i),  dirstr)
               call kwhistp2( rtspec1(j,i),  fno)
            enddo
         enddo
         do j = 1, 3
            do i = 1, hnsites
               call kwhists2( rtspec2(j,i), 0.)
               call kwhistev2(rtspec2(j,i), EventNo)
               call kwhistid2( rtspec2(j,i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")') 
     *              ptcl(j), int( ASDepthList(k)*0.1 )  
               call kseblk(dirstr, "|", nstr)
               call kwhistdir2(rtspec2(j,i),  dirstr)
               call kwhistp2( rtspec2(j,i),  fno)
            enddo
         enddo
      endif


      if(tkretspec ) then
         do j = 1, 3
            do  i = 1, hnsites
               call kwhists3( retspec1(j, i), 0.)
               call kwhistev3(retspec1(j,i), EventNo)
               call kwhistid3( retspec1(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4,"/")' ) 
     *           ptcl(j), int( ASDepthList(k)*0.1 )
               call kseblk(dirstr, "|", nstr)
               call kwhistdir3( retspec1(j, i), dirstr)
               call kwhistp3( retspec1(j, i), fno)
            enddo
         enddo
         do j = 1, 3
            do  i = 1, hnsites
               call kwhists3( retspec2(j, i), 0.)
               call kwhistev3(retspec2(j,i), EventNo)
               call kwhistid3( retspec2(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")' ) 
     *           ptcl(j), int( ASDepthList(k)*0.1 )
               call kseblk(dirstr, "|", nstr)
               call kwhistdir3(retspec2(j, i), dirstr)
               call kwhistp3( retspec2(j, i), fno)
            enddo
         enddo
      endif


      if( tkrezspec ) then
         do j = 1, 3
            do  i = 1, hnsites
               call kwhists3( rezspec(j, i), 0.)
               call kwhistev3(rezspec(j,i), EventNo)
               call kwhistid3( rezspec(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4,"/")' ) 
     *           ptcl(j), int( ASDepthList(k)*0.1 )
               call kseblk(dirstr, "|", nstr)
               call kwhistdir3(rezspec(j, i),  dirstr)
               call kwhistp3( rezspec(j, i),  fno)
            enddo
         enddo
      endif


      if( tkrzfspec ) then
         do j = 1, 3
            do  i = 1, hnsites
               call kwhists3( rzfspec(j, i), 0.)
               call kwhistev3(rzfspec(j,i), EventNo)
               call kwhistid3( rzfspec(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4,"/")' ) 
     *           ptcl(j), int( ASDepthList(k)*0.1 )
               call kseblk(dirstr, "|", nstr)
               call kwhistdir3(rzfspec(j, i),  dirstr)
               call kwhistp3( rzfspec(j, i),  fno)
            enddo
         enddo
      endif


      if( tkrefspec ) then
         do j = 1, 3
            do  i = 1, hnsites
               call kwhists3( refspec(j, i), 0.)
               call kwhistev3(refspec(j,i), EventNo)
               call kwhistid3( refspec(j, i), evid(i))
               k=histdep(i)
               write(dirstr,'(a,"/d",i4, "/")' ) 
     *           ptcl(j), int( ASDepthList(k)*0.1 )
               call kseblk(dirstr, "|", nstr)
               call kwhistdir3( refspec(j, i),  dirstr)
               call kwhistp3( refspec(j, i),  fno)
            enddo
         enddo
      endif
      end
