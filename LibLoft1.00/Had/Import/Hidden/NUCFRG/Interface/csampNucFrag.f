      subroutine csampNucFrag(pj, ita, itz, Nsp, a, np)
!     Form fragments from Nsp spectator nucleons which
!     are left in the collision of  pj and a target of
!     (A,Z)=(ita, itz). If Pj's A is small, not stable
!     in cNucFragXsec
      use modNUCFRAG
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      type(ptcl):: pj  ! input projectile heavy ion
         ! momentum is  not used. only code, subcode,charge
         ! mass and Et are used.

      integer,intent(in):: ita, itz  ! target A, Z
      integer,intent(in):: Nsp   ! # of  spectator nucleons
      type(ptcl)::a(*)  ! output. Nsp nucleons are
                        ! classified into various fragments
                        ! and stored here.  P is given
                        ! assuming  pj is moving z-direction.
      integer,intent(out):: np   ! # of fragmnts in a.
!                 
!  =========================

      integer:: pA, pZ
      real(8):: KEpn
      !      if fragment production Xsec is smaller  than 
      !      this (relative to the total fragment xsec),
      !      neglect it
      real(8),parameter:: minRelXsec=1.0d-4
      
      real(8):: sumNucFrag, sumEmFrag, sumAll
      integer:: nfrag
      integer:: i
      integer:: Zav, dZ,  nnc
      integer:: Zspec, Znuc,  Aspec
      integer:: chg,  Ai
      real(8):: dummyp(3)

      real(8),allocatable::cumSig(:)
      real(8):: u, v
      type(ptcl):: p0
      real(8),parameter:: adjustP = 0.1
      real(8),parameter:: adjustHe = 0.1

      if( pj%code /= kgnuc ) then
         call cerrorMsg('pj is non heavy in csmapNucFrag',0)
      endif
      pA = pj%subcode
      pZ = pj%charge
      KEpn =( pj%fm%p(4) - pj%mass )/pA
!!!!!!!!!!!!
!      write(0,*) ' pA=',pA, ' pZ=',pZ, ' Nsp=',Nsp
!      write(0,*) ' KEpn=',KEpn
!      write(0,*) 'ita, itz =',ita, itz
!!!!!!!/////
      if( KEpn == 0. ) then
         ! target fragment . if 0, abbg shows 
         ! odd result so we put dummy energy
         KEpn = 100E-3
      endif
!         use abbg from cpc.
      call cNucFragXsec(pA, pZ, KEpn, ita, itz) 
!              result is:
!         KOUNT2: # of data in IZZ2,.. SIG3 below
!         IZZ2: Z of fragment
!         IPP2: A of fragment
!         SIG2: fragmentation cross-section (mb)
!         SIG3: photon-induced fragmentation cross-section (mb)
      sumNucFrag = 0.
      sumEmFrag = 0.
      do i = 1, KOUNT2
!!           probably now IsNaN is not needed 
!!          ( at least @JAXA, not usable) so comment out
!!           2014:Mar.6 KK
!!         if( IsNaN(SIG2(i)) ) cycle
!!         if( IsNan(SIG3(i)) ) cycle
         if( SIG2(i) <= 0. .or. SIG3(i) <  0.) cycle
         if( IPP2(i) == 1  .and. IZZ2(i) == 1 ) then
            SIG2(i) = SIG2(i)*adjustP
         elseif( IPP2(i) == 4 .and. IZZ2(i)==2 ) then
            SIG2(i) = SIG2(i)*adjustHe
         endif
         sumNucFrag = sumNucFrag + SIG2(i)
         sumEmFrag =  sumEmFrag + SIG3(i)
      enddo
      sumAll = sumNucFrag + sumEmFrag
!!!!!!!!!!!
!      write(0,*) ' KOUNT2=',KOUNT2
!      write(0,*) ' sumNucFrag =', sumNucFrag 
!      write(0,*) ' sumEmFrag =', sumEmFrag
!      write(0,*) ' sumAll=', sumAll
!!!!!!!!!

!          neglect  fragment of very small XS
      nfrag = 0.
      do i = 1, KOUNT2
!!           see comment for another IsNaN
!!         if( IsNaN(SIG2(i)) ) cycle
!!         if( IsNan(SIG3(i)) ) cycle
         if( SIG2(i) <=  0. .or.  SIG3(i) <  0.) cycle
         if( (SIG2(i) + SIG3(i))/ sumAll  > minRelXsec ) then
            if( IZZ2(i) >0 .and. IPP2(i) > 0 ) then  ! for safety
               nfrag = nfrag + 1
               IZZ2(nfrag) = IZZ2(i)
               IPP2(nfrag) = IPP2(i) 
               SIG2(nfrag) = SIG2(i)
               SIG3(nfrag) = SIG3(i)
            endif
         endif
      enddo

      if(nfrag == 0 ) then
         write(0,*)
     *    'Msg from csampNucFrag: nfrag=0 kount2=',KOUNT2,
     * ' sumAll=',sumAll
         write(0,*) ' pj: A, Z',pj%subcode, pj%charge
         write(0,*) ' pj: KE/n=', (pj%fm%p(4)-pj%mass)/pj%subcode
         write(0,*) ' # of spectator n=', Nsp
         write(0,*) ' target A,Z=',ita, itz
         if( KOUNT2 > 0 ) then
            write(0,*) 'list of sig2,sig3,izz2,ipp2'         
            do i = 1, KOUNT2
               write(0,*) SIG2(i), SIG3(i), IZZ2(i), IPP2(i)
            enddo
         endif
         np = 0
         return !!!!!!!!!!!!!!!!!!
      endif

      sumNucFrag =  sum( SIG2(1:nfrag) )
      sumEmFrag =  sum( SIG3(1:nfrag) )
!!!!!!!!!!
!      write(0,*) ' nfrag=', nfrag
!      write(0,*) ' sumNucFrag=', sumNucFrag
!      write(0,*) ' sumEmFrag=', sumEmFrag
!!!!!!!!!!!!


      allocate(cumSig(nfrag))
!         make cumulative XS table
      cumSig(nfrag) = SIG2(nfrag) + SIG3(nfrag)
      do i = nfrag-1, 1, -1
         cumSig(i) = cumSig(i+1) + SIG2(i) + SIG3(i)
      enddo
!         normalize
      cumSig(:) = cumSig(:)/cumSig(1)
!!!!!!!!!!
!      write(0,*) ' cumSig '
!      do i = nfrag, 1, -1
!         write(0,*) i, cumSig(i), IZZ2(i), IPP2(i)
!      enddo
!!!!!!!!!!!
!           sample fragment; clear counter
      np = 0
      Zspec = 0
      Znuc = 0
      Aspec = Nsp
! 

      do while ( Aspec > 0 )
         call rndc(u)
         do i = nfrag, 1, -1
            if(u <= cumSig(i) ) then
               if(IPP2(i) <= Aspec ) then
                  ! accept this fragment
                  np = np + 1
                  if( IPP2(i) == 1 ) then
                     call cmkptc(knuc,  -1, IZZ2(i), a(np))
                     Znuc = Znuc + 1
                  else
                     if( IZZ2(i) == 2 ) then
                        call rndc(v)
                        if(v < 0.4) then
                           ! split He into 2 deuterons
                           call cmkptc(kgnuc, 2, 1, a(np))
                           np = np + 1
                           call cmkptc(kgnuc, 2, 1, a(np))
                        else
                           call cmkptc(kgnuc, IPP2(i),  IZZ2(i), a(np))
                        endif
                     else
                        call cmkptc(kgnuc, IPP2(i),  IZZ2(i), a(np))
                     endif
                  endif
                  Zspec = Zspec + IZZ2(i)
                  Aspec = Aspec - IPP2(i)
                  if(Aspec == 0 ) exit
                  if(Aspec == 1 ) then
                     np = np + 1
                     call cmkptc(knuc,  -1, 1, a(np))
                     Zspec = Zspec + 1
                     Znuc = Znuc + 1
                     Aspec = 0
                  endif
                  exit
               endif
            endif
         enddo
      enddo
!           assume charge in spectator
      Zav = dble(pj%charge)/pj%subcode *Nsp + 0.5
      if( Zspec > Zav ) then
         dZ = Zspec - Zav
         if( Znuc > dZ ) then
               ! change dZ p to dZ n
            nnc = 0

            do i = 1, np
               if( a(i)%code == knuc ) then
                  nnc = nnc + 1
                  a(i)%charge = 0
                  if(nnc == dZ ) exit
               endif
            enddo
         endif
      endif
!!!!!!!!!!!!
!      do i = 1, np
!         write(0,*)i, a(i)%code, a(i)%subcode, a(i)%charge
!      enddo
!!!!!!

!           give E,P to fragment  and p,n at spectator frame
      do i = 1, np
         chg = a(i)%charge
         if( a(i)%code == knuc ) then
            call csampFermiMom(Nsp, Zav,  chg, a(i)%fm%p)
         else
            Ai = a(i)%subcode
            call csampHFMom(Nsp, Zav, Ai, chg, a(i)%fm%p, dummyp)
         endif
         a(i)%fm%p(4) = 
     *  sqrt( dot_product( a(i)%fm%p(1:3), a(i)%fm%p(1:3)) +
     *       a(i)%mass**2)
      enddo
      p0 = pj
      p0%fm%p(3) = sqrt(pj%fm%p(4)**2 - pj%mass**2)
      p0%fm%p(1:2) = 0.
      do i = 1, np
         call cibst1(i, p0, a(i), a(i))
      enddo
      deallocate( cumSig )
      end subroutine csampNucFrag
