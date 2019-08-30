!       sigma tot, inel
#include "ZcosmosBD.h"
      module modgetXS
      implicit none
      integer,parameter::maxel=15   ! max # of elements
      integer::nel  !  actual @ of elements
      real(8):: A(maxel)
      real(8):: portion(maxel)
!      integer,parameter:: nel=3   ! # of elements
!      real(8):: A(nel)=(/14.0, 16.0, 40.0/)   ! mass # (integer value is used)
!      real(8):: portion(nel)=(/78.1, 20.95, 0.94/)    ! relative # of A's in unit volume

      integer,parameter:: npj=3  ! p, pi, K
      real(8),parameter:: dE=0.0250   ! 0.01 log 10 step
      integer,parameter:: nebin=400  ! # of  Energy bins ( to 10^20 eV) log10
!                                     
      real(8),parameter:: p1= 0.2    ! from p1 GeV/c


      contains
      subroutine computeXS(pj, xsa)
      implicit none
#include "Zptcl.h"
      type(ptcl):: pj  !  projectile . E is fixed here
      real(8),intent(out):: xsa(nebin)
      
      integer:: i, j
      real(8):: p, sum, xs, E, xsI
      real(8):: TA, Z
      p = p1
      do i = 1, nebin
         E = sqrt(p**2 + pj.mass**2)
         pj.fm.p(4) = E
         sum = 0.
         do  j = 1, nel
            TA = A(j)
            Z = max(TA/2, 1.0)
            call cinelx(pj, TA, Z, xsI)
            sum = sum + xsI*portion(j)
         enddo
         xsa(i) =  sum
         p = p*10.0d0**dE
      enddo
      end   subroutine computeXS
      end       module modgetXS
      program main
      use  modgetXS
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zevhnp.h"

      type(ptcl):: pj(3)

      integer:: i, j
      real(8):: xs( nebin,  npj )

      real(8)::  E, temp, s, roots, p
      integer::nela(1)
      
      write(0,*)'current  SxAbySxpOpt',SxAbySxpOpt
      write(0,*)'Enter(1~4) if want to change'
      read(*,*) SxAbySxpOpt
      write(0,*)'Now SxAbySxpOpt',SxAbySxpOpt
      write(0,*)
     * "Enter target A (or A's for compound target) with",
     * "  / at last (<=15)"
      A(:) = 0.
      read(*,*) A(:)
      nela = minloc(A(:))
      nel = nela(1) -1
      if( nel == 0 ) stop
      if( nel == 1 ) then
         portion(nel) = 1.
      else
         write(0,*) "Enter relative portion of A's"
         read(*,*) portion(1:nel)
      endif

      temp = sum(portion(1:nel))
      portion(:) = portion(:)/temp


      ! make proton
      call cmkptc(knuc, -1, 1, pj(1))
      call computeXS(pj(1), xs(1, 1))  ! p
      ! make pi
      call cmkptc(kpion, -1, 1, pj(2))
      call computeXS(pj(2), xs(1, 2))  ! pi
      ! make kaon
      call cmkptc(kkaon, -1, 1, pj(3))
      call computeXS(pj(3), xs(1, 3))  ! K
!         output

      do j = 1, npj 
         p= p1
         do i = 1, nebin
            E = sqrt(p**2 + pj(j).mass**2)
            s = pj(1).mass**2 + pj(j).mass**2 + 2.*E*pj(1).mass
            roots = sqrt(s)
            write(*,'(a,i1, 1p, 5g13.4 )')
     *        'pj ',j, E, xs(i, j), p, s, roots
            p = p*10.d0** dE
         enddo
      enddo
      end program main

