      module modgetXsRatio
      implicit none
      integer,parameter:: nel=43
      real(8):: A(nel)=(/1.,2.,4.,5.,6.,8.,10.,12.,14.,16.,18.,20.,
     * 22.,24.,26.,28.,30.,32.,36.,40.,
     * 44.,48.,52.,56.,60.,64.,68.,72.,76.,80.,90.,100.,
     * 110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210./)

      integer,parameter:: npj=3  ! p, pi, K
      integer,parameter:: nebin=110  ! # of  Energy bins (from  10^10 eV to 10^21 eV) log10
!                                     ! step 0.1
      real(8),parameter:: E1= 10.    ! from 10GeV
      integer,parameter:: Eout =   70.  ! E> Eout will be printed
      real(8),parameter:: dE=0.1d0   ! 0.1 log 10 step

      contains
      subroutine computeXS(pj)
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj  !  projectile . E is fixed here

      type(ptcl)::tg 
      integer:: i, j
      real(8):: E,  xsp, xs
      integer TA
      E = E1

      do i = 1,  nebin
         if(E > Eout) then
            pj%fm%p(4) = E
            TA= A(1)            !make proton target
            call cmkptc(knuc, -1, 1, tg) 

            call  cgetXsInterface2( pj, tg,   xsp )
            do  j = 2, nel
               TA = A(j)
               call cmkptc(kgnuc, TA, max(TA/2,1), tg)
               call  cgetXsInterface2( pj, tg,   xs )
               write(*,'(i3, 1p,g14.4, 0p, i4, f9.1, f9.1, f8.4)')
     *            pj%code, pj%fm%p(4),  TA, xsp, xs,  xs/xsp
            enddo
         endif
         E = E*10.0d0**dE
      enddo
      end   subroutine computeXS
      end       module  modgetXsRatio
      program main
      use  modgetXsRatio
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl):: pj(3)

      integer:: i, j


      call cgetXsIni
      write(0,'(a)')
     *  '# pj   E(GeV)  Tartget A s(xp)(mb)  s(xA) s(xA)/s(xp) '
      ! make proton
      call cmkptc(knuc, -1, 1, pj(1))
      call computeXS(pj(1))
      ! make pi
      call cmkptc(kpion, -1, 1, pj(2))
      call computeXS(pj(2))
      ! make kaon
      call cmkptc(kkaon, -1, 1, pj(3))
      call computeXS(pj(3))
      
      end program main
