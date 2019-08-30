      subroutine cfrevent(kf, charge, ia, iz, roots,  a, ntp)
!
      implicit none

#include  "Zptcl.h"
#include "Zkfcode.h"
#include "Zmanagerp.h"
#include "Zevhnp.h"
! &&&&&&&&&&
#include "Zcode.h"
#include "Zmass.h"
!  &&&&&&&&&&

!     *******************
      include "Zlujets.h"
!     *******************

!
      integer kf    ! input. projectile kf code
      integer charge ! input. projectile charge
      integer ia     ! input. nucleon number of the target
      integer iz     ! input. charge no. of target
      real*8 roots  ! input. cms energy (GeV)
      type (ptcl):: a(*)  !  output. produced ptcls.  in cms.
      integer ntp   ! number of produced ptcls
!
      logical first/.true./
      save first

      integer code, subcode, chg, i

      integer ipt,  nnuc, nprot, kcd
      real*8  ro1, exma
      character*4 pacd 
      common/frcodes/ipt(2),pacd(27),nnuc(27),nprot(27),kcd(27)
     >           ,ro1(27,2),exma(9,2)

      integer mstp, msti
      real*8 parp, pari
      common/pyparsC/mstp(200),parp(200),msti(200),pari(200)
!                   
      integer ksz1
      parameter (ksz1=20)
      real*8 vfr
      integer kfr
      common/frpara1/kfr(ksz1),vfr(ksz1)

      integer mstu, mstj
      real*8 paru, parj
      common/ludat1/mstu(200),paru(200),mstj(200),parj(200)
      integer mdcy, mdme, kfdp
      real*8 brat
      common/ludat3/mdcy(500,3),mdme(2000,2),brat(2000),kfdp(2000,5)
      real*8 crmsnc, cr0Woodsaxon, ccWoodsaxon

      character*4 pjid, tgtid


      save /ludat1/,/lujets/,/ludat3/

      integer lucomp, ind
!////////////
!      integer seed(2)
!////////////

      if(first) then
!           lujets size
         mstu(4) = kszj
!           no output from Fritiof
         mstu(12) = 0
!         compress lujets
         kfr(13) = 1
         kfr(11) = LundPara(4)  !  write when parameter change to frevent
         kfr(12) = LundPara(2)  ! opal or delphi
         mstp(127) = LundPara(3)  ! suppress pythia message
         mstp(122) = LundPara(3)  ! //
         mstu(11) = ErrorOut
!
!         allow the decay of lambda
         mdcy(lucomp(3122),1) = 1
!          suppress k_s0:  decay
         mdcy(lucomp(111), 1) = 0
!          forbid pi0 decay
         mdcy(lucomp(310),1) = 0
         first = .false.
      endif
      if(roots .gt. 4.5d5) then
         mdcy(lucomp(221), 1)= 0 ! forbid eta decay
      else
         mdcy(lucomp(221), 1)= 1 ! allow eta decay
      endif
!         use pythia multiple hard rps.(2) or simple fritiof hard scat.
!         (=1) (this may be changed in Fritiof to 0 at low energy,
!         so don't put it in the above init part. 
!        
      kfr(7) = LundPara(1)
!          set projectile and target info.
      call ckf2pacd(kf, ind) ! kf to pacd array index. 0 is no corres.
      if(ind .eq. 0) then
         pjid ='new1'
         kcd(1) = kf
         nprot(1) = charge
         ipt(1) = 1  ! may not be needed
         if(kf .eq. kfk0l .or. kf .eq. kfk0s) then
            exma(1,1) = exma(5,1)
            exma(1,2) = exma(5,2)
!         elseif(kf .eq. kfnbar) then
         elseif(kf .eq. -kfneutron) then
            exma(1,1) = exma(9,1)
            exma(1,2) = exma(9,2)
         elseif(kf .eq. kfpi0) then
            exma(1,1) = exma(3,1)
            exma(1,2) = exma(3,2)
         else
            write(*, *) 'strange kf code=', kf, ' to  cfrevent'
            stop 9999
         endif

      else
         pjid = pacd(ind)
      endif

!             set target info.
      call caz2pacd(ia, iz, ind)     ! target index in pacd array.
!     
      if(ind .eq. 0) then
         tgtid ='new2'
         nnuc(2) = ia
         nprot(2) = iz
         ipt(2) = 2   ! may not be needed
         if(ia .le. 16) then
            ro1(2,1) = crmsnc(iz)
         else
            ro1(2,1) = cr0Woodsaxon(iz)
            ro1(2,2) = ccWoodsaxon(iz)
         endif
      else
         tgtid = pacd(ind)
      endif
!    ---------------
      n = 2
      do while( n .eq. 2)       ! avoid elestic events
!//////////////////
!         if(roots .gt. 1.5e2) then
!            write(*,*) ' ====== roots=',roots,
!     *            pjid, tgtid
!            call rnd1s(seed)
!            write(*,*) ' seed=', seed
!         endif   
!////////////////
         call frevent('cms', pjid, tgtid, roots)
      enddo

!       data in /lujets/ is moved to 'a' with Cosmos code
      ntp = 0
      do i = 1, n
         call ckf2cos(k(i, 2), code, subcode, chg)
         if(code .gt. 0) then
            ntp = ntp + 1
            call cmkptc(code, subcode, chg, a(ntp))
            a(ntp)%fm%p(1) = p(i, 1)
            a(ntp)%fm%p(2) = p(i, 2)
            a(ntp)%fm%p(3) = p(i, 3)
            a(ntp)%fm%p(4) = p(i, 4)
         else
!   &&&&&&&&&&&&&&&&&&&&
!             though rare, strange ptcl comes out with charge << -1
!             and mass = 0. . neglect  such one
!               totally neglect
!!!            if( p(i,5) .le. 0.) goto 10
!!!            subcode = p(i,5)/masp +  0.5
!!!            chg = abs( k(i,2) ) - 10000
!!!            if(chg  .le. 0 ) goto 10
!!!            code = kgnuc
!!!c            &&&&&&&&&&&&&& 
!!!            if(subcode .eq. 0) then
!!!c              sometimes Fritiof forgets to set mass.
!!!c              assume it is 2.2*charge
!!!               subcode = chg*2.2
!!!            endif
!!!            ntp = ntp + 1
!!!            call cmkptc(code, subcode, chg, a(ntp))
!!!            a(ntp).fm.p(1) = p(i, 1)
!!!            a(ntp).fm.p(2) = p(i, 2)
!!!            a(ntp).fm.p(3) = p(i, 3)
!!!            a(ntp).fm.p(4) = p(i, 4)
!!!            a(ntp).mass = p(i, 5)
!!!c &&&&&&&&&&&&&&&&&&
         endif
 10      continue
      enddo
      end
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ckf2pacd(kf, ind)
      implicit none
#include "Zkfcode.h"
!      get index for pacd array from kf code
      integer kf ! input kf code.
      integer ind ! output.  kf code corresponds to ind-th position
!                        in pacd array in Fritiof. 0 is given if
!                        no correspondence.
!       For heavy partilce, 0 is returned. (only other hadrons are 
!       treated as projectile in Cosmos.  heavy particle is
!       not directly treated by Lund.

      if(kf .eq. kfpion) then
         ind = 3
      elseif(kf .eq. -kfpion) then
         ind = 4
      elseif(kf .eq. kfkaon) then
         ind = 5
      elseif(kf .eq. -kfkaon) then
         ind = 6
      elseif(kf .eq. kfneutron) then
         ind = 7
      elseif(kf .eq. kfproton) then
         ind = 8
      elseif(kf .eq. -kfproton) then
         ind = 9
      else
         ind = 0
      endif
      end
      subroutine caz2pacd(ia, iz, ind)
      implicit none
!        given nucleus is converted pacd array index in Fritiof
      integer ia  ! input. nucleon numbr
      integer iz  ! input. proton number
      integer ind ! output. index for pacd arrary.


      integer ipt,  nnuc, nprot, kcd
      real*8  ro1, exma
      character*4 pacd 
      common/frcodes/ipt(2),pacd(27),nnuc(27),nprot(27),kcd(27)
     >           ,ro1(27,2),exma(9,2)

!
      integer i
          

      if(ia .eq. 1  .and. iz .eq. 1) then
            ind = 8  ! proton;  n, p_bar cannot be target in Fritiof
      else
         do i =10, 27
            if(nprot(i) .eq. iz) then
               ind = i
               goto 100
            endif
            ind = 0
         enddo
 100     continue
      endif
      end
!             test cr0Woodsaxon
!	integer z
!       real*8 cr0Woodsaxon
!
!	do z= 8, 100, 2
!	    write(*, *) z, sngl(cr0Woodsaxon(z))
!	enddo
!	end

         real*8 function cr0Woodsaxon(z)
         implicit none
!            	
        integer z       ! input atomic number > 8
         integer i
         real*8 ans

         real*8 r0(6) ! radius  parameters for woodsaxon density
!                     for Z> 8 
          data ( r0(i), i=  1,   6)/
     1   0.85744508    ,  0.17366389E-01, -0.54510422E-03,
     2   0.91348304E-05, -0.76877989E-07,  0.25461102E-09 
     * /   

         

         
         ans = 0.
         do i = 6, 1, -1
            ans = ans*z + r0(i)
        enddo	
         cr0Woodsaxon = ans
         end
!             test ccWoodsaxon
!	integer z
!       real*8 ccWoodsaxon
!
!	do z= 8, 100, 2
!	    write(*, *) z, sngl(ccWoodsaxon(z))
!	enddo
!	end

         real*8 function ccWoodsaxon(z)
!              compute edge parameter of WoodSaxon density.
         implicit none
!            	
        integer z       ! input atomic number > 8
         integer i
         real*8 ans

         real*8 c(10) ! edge parameters for woodsaxon density
!                     for Z> 8 
          data ( c(i), i=  1,   10)/
     1   -1.2235859    ,  0.45177680    , -0.50751763E-01,
     2         0.31677493E-02, -0.12043819E-03,  0.28859657E-05,
     3  -0.43630795E-07,  0.40272744E-09, -0.20689960E-11,
     4         0.45292573E-14                                                 
     * /   

         if(z .lt. 13) then
            ccWoodsaxon = 0.47
         else
            ans = 0.
            do i = 10, 1, -1
               ans = ans*z + c(i)
            enddo	
            ccWoodsaxon = ans
         endif
         end
      real*8 function crmsnc(z)
!       root mean square radius 
      implicit none
      integer z ! input. nuclear charge
!     
      real*8 rms(8)  ! for z=1,2,3,4,5,6,7,8 (d,He,Li,Be,B,C,N,O)  
      data rms/
     1     2.095, 1.74,  2.16, 2.519, 2.37, 2.466, 2.52, 2.724/

      crmsnc= rms(z)
      end
