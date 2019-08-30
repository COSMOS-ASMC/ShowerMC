#include "BlockData/cblkGene.h"
      module modgetDPMxs
      implicit none
      integer,parameter:: nel=1   ! # of elements
      real(8):: A(nel)=(/1.0/)
      real(8):: portion(nel)=(/1./)

!      integer,parameter:: nel=3   ! # of elements
!      real(8):: A(nel)=(/14.0, 16.0, 40.0/)   ! mass # (integer value is used)
!      real(8):: portion(nel)=(/78.1, 20.95, 0.94/)    ! relative # of A's in unit volume
      integer,parameter:: npj=3  ! p, pi, K
      integer,parameter:: nebin=110  ! # of  Energy bins (from  10^10 eV to 10^21 eV) log10
!                                     ! step 0.1
      real(8),parameter:: E1= 10.    ! from 10GeV
      real(8),parameter:: dE=0.1d0   ! 0.1 log 10 step

      contains
      subroutine computeXS(pj, xsa)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"

      type(ptcl):: pj  !  projectile . E is fixed here
      real(8),intent(out):: xsa(nebin)

      integer:: i, j
      real(8):: E, sum, xs, mfp
      integer TA

      TrackBefMove.p = pj
      E = E1
      do i = 1, nebin
         TrackBefMove.p.fm.p(4)= E
         sum = 0.
         do  j = 1, nel
            TargetMassN = A(j)
            call cmfpdpmjet3(xs, mfp)
            sum = sum + xs*portion(j)
         enddo
         xsa(i) =  sum
         E = E*10.0d0**dE
      enddo
      end   subroutine computeXS
      end       module modgetDPMxs
      program main
      use  modgetDPMxs
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj

      integer:: i, j
      real(8):: xs( nebin,  npj )

      real(8)::  E
      portion(:) = portion(:)/sum(portion(1:nel))

      ! make proton
      call cmkptc(knuc, -1, 1, pj)
      call computeXS(pj, xs(1, 1))  ! p
      ! make pi
      call cmkptc(kpion, -1, 1, pj)
      call computeXS(pj, xs(1, 2))  ! pi
      ! make kaon
      call cmkptc(kkaon, -1, 1, pj)
      call computeXS(pj, xs(1, 3))  ! K
!         output
      E= E1
      do i = 1, nebin
         write(*,'(1p, 4g13.4 )')  E, (xs(i, j), j= 1, npj)
         E = E*10.d0** dE
      enddo
      end program main




