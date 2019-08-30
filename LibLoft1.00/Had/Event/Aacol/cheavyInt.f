!     *****************************************************************
!     *                                                               *
!     *  cheavyInt: treats heavy particle interactions
!     *                                                               *
!     *****************************************************************
!
!  -- process --
!             1) samples fragment ptcls using fragmentation parameters
!                and determines the no. of interacting nucleons.
!             2) gives break-up angle for fragments other than interacting
!                nucleons
!             3) for interacting nucleons, calls  chAcol to make
!                multiple production.
!
        subroutine cheavyInt(pj, ia, iz, xs, a, ntp)
! #if      !defined (MacIFC)
!    same  note as in chAcol.f     
!        use modXsecMedia
! #endif

!        use modColInfo
        implicit none

#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zcoord.h"
#include  "Zheavyv.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"


!///////////
!        logical deb
!        common/cdebug/ deb
!//////////
!
        external cblkHeavy
!
!
!          ia: mass no. of the target
!          iz: charge of the target

        integer ia, iz, ntp
        real(8),intent(in):: xs  ! x-section on A= ia mb
        type(ptcl):: pj, a(*) 
!
!        type(coord):: dir
!        real*8  dummylen
!
!        to store heavy fragment            nucleon
        type(ptcl):: frga(maxHeavyMassN), nuc(maxHeavyMassN)

        type(ptcl):: fragA(maxHeavyMassN), intNucA(maxHeavyMassN),
     *    nonIntNucA(maxHeavyMassN)
!

        integer noOfFragments, noOfNuc, noOfInteNuc
        save  noOfFragments, noOfNuc, noOfInteNuc

        integer noOfFrag, noOfIntN, noOfNonIntN

        integer  i, n, j
!///////////////
!        integer cumev, cev
!////////////
!         /////
!         real*8 sumerg

!         call ctestOnShell(' heavy bef frag', pj)
!         sumerg = 0.
!         ////

!     
! MacIFC  ld error  ; also for TargetXs here and chAcol.f)
!  (No error when complile ephook.f but for
!   testCnf1.mk error occurs)
!////        integer:: dummy
!///         dummy = TargetNucleonNo
!Undefined symbols for architecture x86_64:
!  "_modxsecmedia_mp_targetnucleonno_", referenced from:
!      _cheavyint_ in libcosmos.a(cheavyInt.o)
!  "_modxsecmedia_mp_targetxs_", referenced from:
!      _cheavyint_ in libcosmos.a(cheavyInt.o)
!      _chacolx_ in libcosmos.a(chAcol.o)
!ld: symbol(s) not found for architecture x86_64
!c//////////////////
#if defined (MacIFC)
!     #include "Zworkaround.h"
#endif


        
      ntp=0

!
!            ** fragmentation with some nucleon interaction **
!                 sample and set fragmentation ptcls
!      noOfFragments:  heavy fragments
!      noOfNuc:  all nucleons
!      noOfInteNuc: interacting nucleons.
!
!                                          
!        if(IntModel .eq. 'int2' .and.
!     *      ( pj.code .eq. ktriton .or. pj.code .eq. kdeut) ) then

        if( ActiveMdl .eq. 'dpmjet3') then
           if((pj%fm%p(4) - pj%mass)/pj%subcode .gt. 5.1) then
              call cdpmjet(pj, ia, iz, a, ntp)
           else
!               rescue
              call cjamEvent(pj, ia, iz, xs, a, ntp)
           endif
        elseif( ActiveMdl .eq. 'jam') then
           call cjamEvent(pj, ia, iz, xs, a, ntp)
        elseif( ActiveMdl .eq. 'phits') then
           call cphits(pj, ia, iz, xs, a, ntp)
        elseif(ActiveMdl .eq. 'qgsjet2') then
           call cQGSjet(pj, ia, iz, a, ntp)
        elseif(ActiveMdl .eq. 'epos') then
           call ceposGenOneEvent(pj, ia, iz, a, ntp)
        elseif(ActiveMdl .eq. 'gheisha' .and.
     *        pj%code .eq. kgnuc  .and. (pj%subcode .eq. 2 .or.
     *        pj%subcode .eq. 3) ) then
!                  Gheisha should be called directly
           call chAGheisha(pj, ia, iz, a, ntp)
           noOfFragments = 0
           noOfNuc = 0
           noOfInteNuc = 0
        else
!              resque
           call cjamEvent(pj, ia, iz, xs, a, ntp)           
!!           call csampFragments(pj, ia, frga, noOfFragments,
!!     *       nuc, noOfNuc, noOfInteNuc)
!!           if(SkipPtclGen .eq. 0) then
!!c                 interaction of interating nuc.
!!              do  i = 1, noOfInteNuc
!!                 j=i + noOfNuc - noOfInteNuc
!!                 call chAcol(nuc(j), ia, iz, a(ntp+1), n)
!!                 ntp=ntp+n
!!              enddo
!!c               all business so far is done in the frame where
!!c               z axis is the incident
!!c               make rotation
!!              call crot3mom( pj, a, ntp )
!!           endif
!!           call crot3mom( pj, nuc, noOfNuc) ! for all nucleons. see entry below
!!
!!           call crot3mom( pj, frga, noOfFragments)
!!
!!
!!c               move non-interacting nucleon and fragments
!!           do  i = 1, noOfNuc - noOfInteNuc 
!!              a(ntp+1) = nuc(i)
!!              ntp = ntp +1
!!           enddo
!!           do  i = 1, noOfFragments
!!              a(ntp + 1) = frga(i)
!!              ntp = ntp + 1
!!           enddo
        endif
        return
!     ****************** inquire the fragments at heavy interation
      entry cqHvyIntF(fragA, noOfFrag)
!     ******************* 
!        note if MovedTrack.p.code is not heavy, this gives wrong result

      noOfFrag = noOfFragments
!         move fragments
      do i = 1, noOfFrag
         fragA(i) = frga(i)
      enddo
      return
!     ******************** inquire the interacting nuc.
      entry cqHvyIntIN(intNucA, noOfIntN)
!     *******************
      noOfIntN = noOfInteNuc
!        move  interacting nucleons
      do i = 1, noOfIntN
         intNucA(i) = nuc(i+noOfNuc - noOfInteNuc)
      enddo
      return
!     ******************inquire non interacting  nuc.
      entry cqHvyIntNIN(nonIntNucA,  noOfNonIntN)
!     ******************
      noOfNonIntN = noOfNuc - noOfInteNuc
!        move non-interacting nucleons
      do i = 1, noOfNonIntN
         nonIntNucA(i) = nuc(i)
      enddo
      end
!     *****************************************************************
!     *                                                               *
!     *  csampFragments: sample and set fragmentation ptcls           *
!     *                                                               *
!     *****************************************************************
!
!      /usage/  call csampFragments(pj, ia, fra, noOfFragments, 
!    *    nuc,  noOfNuc, noOfInteNuc)
!
!       1)  samples fragment ptcl one by one by referring cfptbl
!           until sum of fragment mass no. exceeds incident mass no.
!           - 3 so that no heavy ptcl can emerg more.  then the 
!           remaining ones
!           are assigned to nucleons, if any.  If the sum exceeds
!           incident mass no. during sampling, retrial is made from the
!           first.   For each sampled fragment, energy is given pro-
!           potionally to its mass.   charge, mass and kind are also
!           given according to cosmos convention.
!
!        2) fra is to store heavy fragments
!           nuc is to store all nucleons
!        3) charge of the nucleons is reset because process 1) assigns
!           only proton
!        4) samples no. of interacting nucleons
!        5) samples pt of fragments other than interacting nucleons
!           and convert it to ptx, pty
!
!      noOfFragments:  # of heavy fragments stored in fra
!      noOfNuc:       # of nucleons  stored in nuc.
!      noOfInteNuc:  # of interacting nucleons among them
!
!
        subroutine csampFragments(pj, ia,
     *       fra,  noOfFragments, nuc, noOfNuc, noOfInteNuc)
        implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zheavyp.h"
#include  "Zheavyc.h"


        integer noOfFragments, noOfNuc, noOfInteNuc, ia
!         ia: # of target nucleons
        type(ptcl):: pj, fra(*), nuc(*)
        integer ihg, mno, mx, mn, msumf, ihgf, jcon,
     *          hchg
        real*8  epn, u
        logical first/.true./
        integer i, zfrag
!
      if(first) then
         call cmakeFragTbl
         first = .false.
      endif

      if( (pj%fm%p(4)-pj%mass)/pj%subcode .lt. 0.1) then
         fra(1) = pj
         noOfFragments = 1
         noOfNuc = 0
         noOfInteNuc = 0
         return   !   ********
      endif
!
!         get heavy group index from the incident charge  
      ihg=Charge2heavyG(pj%charge)

!          # of nucleons
!      mno=HeavyG2massN(ihg)
      mno = pj%subcode
!       z 
!      hchg=HeavyG2charge(ihg)

       hchg=pj%charge

!       t.e  energy per nucleon
      epn=pj%fm%p(4)/mno
!         margin for mass no. conservation (for first trial)
      mx = mno
      mn = mno - 3

!
!             retry if final sum of mass exceeds incident one
!      *** until loop*** 
      do while (.true.)
!            sum of mass of fragments
         msumf = 0
!            no. of heavy fragments
         noOfFragments = 0
!            no. of all nucleons 
         noOfNuc =  0
!               repeat until sum of mass becomes >= amn
!         *** until loop*** 
         do while (.true.)
            if( mno .ge. 10  .and. noOfFragments .eq. 0) then
!                   for A >= 10, first sampled fragment must be
!                   non-nucleon to avoid too many nucleons.
                  call rndc(u)
                  u=(1.-CfragmentationTbl(ihg,2))*u +
     *                  CfragmentationTbl(ihg,2)
!                     find CfragmentationTbl >= u to sample fragrment
                  call kfrge(CfragmentationTbl(ihg,1), maxHeavyG,
     *            ihg, u,  ihgf, jcon)
!                            sampled group index
            else
               call rndc(u)
!                    find FragmentaionTbl >= u to sample fragrment
               call kfrge(CfragmentationTbl(ihg,1), maxHeavyG,
     *           ihg, u,  ihgf, jcon)
            endif
!
            if(ihgf .eq. 1)then
!               count nucleon
               noOfNuc = noOfNuc + 1
!               make ptcl;  nucleon charge is reset later
               call cmkptc(knuc, regptcl,
     *          1, nuc(noOfNuc))
               nuc(noOfNuc)%fm%p(4) = epn
            else
!               count fragment
               noOfFragments=noOfFragments+1
!               make ptcl; 
               call cmkptc(HeavyG2code(ihgf), regptcl,
     *          1, fra(noOfFragments))
!                 set energy
               fra(noOfFragments)%fm%p(4) = HeavyG2massN(ihgf) * epn
            endif
!               count mass no.
            msumf=msumf+HeavyG2massN(ihgf)
            if  (msumf .ge. mn)
     *                      goto 50
         enddo
   50    continue
         if(msumf .le. mx) then
!           remaining ptcls should be nucleons
            do i = 1, mno - msumf
               noOfNuc = noOfNuc + 1
               call cmkptc(knuc, regptcl, 1, nuc(noOfNuc) )
!                set energy
               nuc(noOfNuc)%fm%p(4) = epn
            enddo
!             get sum of fragment charge
            zfrag = 0
            do i = 1, noOfFragments
               zfrag = zfrag + fra(i)%charge
            enddo
            if(zfrag .le. pj%charge) goto 100
         endif
      enddo
  100 continue
!          reset nucleon charge
      call cresetNucChg(nuc, noOfNuc, hchg - zfrag)
!
!               sample interacting nucleon no.
      call csampInteNuc(pj, ia,  noOfNuc, noOfInteNuc)

!           sample fragment mom. all frag. 
      call csampFragMom(fra,  noOfFragments)
!           sample non-interacting nuc. mom.
      call csampFragMom(nuc, noOfNuc - noOfInteNuc)
!           set interacting nuc. mom
      do   i=noOfNuc, noOfNuc - noOfInteNuc + 1, -1
         nuc(i)%fm%p(1) = 0.
         nuc(i)%fm%p(2) = 0.
         nuc(i)%fm%p(3) = sqrt(
     *      max( nuc(i)%fm%p(4)**2 - nuc(i)%mass**2, 0.d0) )
      enddo
      end
!     *****************************************************************
!     *                                                               *
!     *  cresetNucChg:  reset charge of nucleons emerging from heavy *
!     *                                                               *
!     *****************************************************************
!
       subroutine cresetNucChg(nuc, noOfNuc, z)
!
       implicit none
!----       include '../../Zptcl.h'
#include  "Zptcl.h"
       integer z, noOfNuc
       type(ptcl):: nuc(noOfNuc)
!
       integer i
!
       do i = 1, z
          nuc(i)%charge = 1
       enddo
       do i = z+1, noOfNuc
          nuc(i)%charge = 0
       enddo
      end
!     *****************************************************************
!     *                                                               *
!     *  csampInteNuc:  sample interacting nucleon number                   *
!     *                                                               *
!     *****************************************************************
!
!        1) gets average no. of interacting nucleons
!        2) gets average no. of nucleons from heavy
!        3) using the ratio of 1)/2) and binomial distribution,
!           samples no. of interacting nucleons.

!
!                           =   =   =   =
!
      subroutine csampInteNuc(pj, ia,  noOfNuc, noOfInteNuc)

      implicit none

#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zheavyp.h"
      
      integer noOfNuc, noOfInteNuc, ia  ! ia: # of target nuc. 
      type(ptcl):: pj  ! projectile heavy
      integer ihg
      real*8  avintn, avnn
!

!                heavy group index
      ihg=Charge2heavyG(pj%charge)
!                get average no. of interacting nucleons
      call caveInteNuc(pj, ia,  avintn)
!                average no. of nucleons in fragments
      avnn = FragmentTbl(ihg, 1)
!       sample interacting nucleon number by binormial distribution
!       with prob. avintn/avnn
      call kbinom( avintn/avnn, noOfNuc, noOfInteNuc)
      end
!     *****************************************************************
!     *                                                               *
!     *  csampFragMom:  sample  px, py, pz of fragments
!     *                                                               *
!     *****************************************************************
!
!      /usage/     call csampFragMom(a, nf)
!
!        pt is sampled by gaussian type distribution and stored in
!        'a'. ptx,pty are also stored.
!
      subroutine csampFragMom( a, nf )
      implicit none

#include  "Zptcl.h"
      integer nf
      type(ptcl):: a(nf)

      integer i, nc
      real*8 pt, p, cs, sn

       do   i=1, nf
!            sample fragment pt
         nc=0
         p=sqrt( a(i)%fm%p(4)**2- a(i)%mass**2 )
!         *** until loop*** 
         do while (.true.)
            call csampFragPt(a(i),  pt)
            nc=nc+1
            if         (pt .lt. p .or. nc .eq. 10)
     *                      goto 10
         enddo
   10    continue
         if(nc .ge. 10) then
             pt=min(1.d-10, p)
         endif
!               set pt and pz
         a(i)%fm%p(3) = sqrt(p**2-pt**2)
         call kcossn(cs, sn)
         a(i)%fm%p(1) = pt*cs
         a(i)%fm%p(2) = pt*sn
       enddo
      end
!     *****************************************************************
!     *                                                               *
!     *  csampFragPt:   sample fragment pt                      
!     *                                                               *
!     *****************************************************************
!
!  -- process --
!        pt is sampled by gaussian type distribution:
!           exp(- x**2) dx**2  where x= 2pt/sqrt(pi)
!
      subroutine csampFragPt(aPtcl, pt)
      implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zheavyp.h"
      real*8  pt      ! sampled pt in GeV/c

      type(ptcl):: aPtcl
!
!
      real*8 avpt, u
!
      if(aPtcl%code .eq. knuc) then
        avpt= PtAvNonInteNuc
      else
        avpt= PtAvFrag
      endif
!
      call rndc(u)
!          note avpt is not <pt> but <pt>/sqrt(pi/2)
      pt  =  avpt *sqrt(- log(u)* 2 )
      end
!     *****************************************************************
!     *                                                               *
!     *  cmakeFragTbl:  makes fragmentation parameter table for 
!     *                 heavy ptcl frgmentation sampling              
!     *                                                               *
!     *****************************************************************
!
!        FragmentaionTbl(i,j) is assumed to have <no. of heavy of group 
!        j>  when heavy of group i fragments.  For each group i, 
!        FragmentTbl is normalized so that the sum of them
!        becomes 1 and then made to be cumulative table such that
!        CfragmentaionTbl(i,j) <= CragmentaionTbl(i,j+).
!
!

      subroutine cmakeFragTbl
      implicit none

#include  "Zcode.h"
#include  "Zheavyp.h"
#include  "Zheavyc.h"
!
      integer i, j
!
!          do below for nucleus group 1 to maxHeavyG
!
       do   i=1, maxHeavyG
!              FragmentaionTbl(i,j) containes 
!              <no. of nucleus of group j> when
!              nucleus of group i fragments  (j<=i)
!                make cumulative table
          CfragmentationTbl(i, 1) = FragmentTbl(i, 1)
          do   j=1,i-1
             CfragmentationTbl(i,j+1)=FragmentTbl(i,j+1)
     *          +CfragmentationTbl(i,j)
          enddo
!                 normalzie
          do  j=1, i  
             CfragmentationTbl(i, j) = CfragmentationTbl(i,j)/
     *          CfragmentationTbl(i,i)
          enddo
       enddo
      end
      subroutine chg2massN(hg, massnum)
      implicit none
#include "Zcode.h"
#include "Zheavyp.h"
      integer massnum, hg
      massnum = HeavyG2massN(hg)
      end
      subroutine chg2charge(hg, charge)
      implicit none
#include "Zcode.h"
#include "Zheavyp.h"
      integer charge, hg
      charge = HeavyG2charge(hg)
      end
      subroutine ccode2hvgrp(code, hg)
      implicit none
#include "Zcode.h"
#include "Zheavyp.h"
      integer code, hg
      if(code .ge. kalfa .and. code .le. kiron) then
         hg = Code2heavyG(code)
      else
         call cerrorMsg(
     *   'ccode2hvgrp should not be used for code # He ~Fe',0)
      endif
      end
      subroutine ccode2mass(code, mass)
      implicit none
#include "Zcode.h"
#include "Zheavyp.h"
#include "Zmass.h"
      integer code, hg
      real*8 mass
      if(code .ge. kalfa .and. code .le. kiron) then
         hg = Code2heavyG(code)
         mass = HeavyG2massN(hg) * (masp+masn)/2
      else
         call cerrorMsg(
     *   'ccode2mass should not be used for code # He ~Fe',0)
      endif
      end

      
