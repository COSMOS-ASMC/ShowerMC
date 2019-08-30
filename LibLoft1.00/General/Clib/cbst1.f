#include "ZsaveStruc.h"
!          testing cbst0, 1
!     include 'cbst0.f'
!     include 'clorep.f'
!     include 'clorez.f'
!     include 'cgetRotMat4.f'
!     include 'cgetlf.f'
!     include 'cpm2e.f'
!     include 'cibst1.f'
!c-----------------------
!      implicit none
!      include '../Zptcl.h'
! 
!       type(ptcl):: p1, p2, po, pox
!       real*8 m
!
!       m=.938
!
!           p1.fm.p(1)=1.5
!           p1.fm.p(2)=-8.5e1
!           p1.fm.p(3)=5.e4
!           p2.fm.p(1)=-1.0
!           p2.fm.p(2)=10.
!           p2.fm.p(3)=3000.
!           p1.mass=m
!           p2.mass=m
!           call cpm2e(p2, p2)
!           call cpm2e(p1, p1)
!           call cbst1(1,  p1, p1, po)
!           call cibst1(1, p1, po, pox)
!           write(*,*) po.fm.p, pox.fm.p
!       end
!       **************************************************************
!       *
!       * cbst1:  boost a partilce into the rest system of another
!       *         particle
!       *
!       **************************************************************
!
! /usage/  call cbst1(init, p1, p2,  po)
!
!        Suppose two particles p1, p2, given in the system K. 
!        This boosts particle 2 into  the rest frame of particle 1.
!        The z axis of p1 is made to be the direction of  p1.
!  init: integer.  Input.  If p1 is the same as the previous call to this
!        subroutine, give a value other than 1.  If p1 is different
!        from the previous call, use 1. 
!    p1: type ptcl. Input. particle 1
!    p2: type ptcl. Input. particle 2
!    po: type ptcl. Output. particle seen at the rest system
!        of particle 1.
!
       subroutine cbst1(init, p1, p2, po)
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
         type(ptcl):: p1, p2, po
         integer init
!
         type(fmom):: g
#ifdef USESAVE
         save g
#endif
!
         if(init .eq. 1) then
!                 get Lorent factor of p1
            call cgetlf(p1, g)
         endif   
!              boost
         call cbst0(init, g, p2,  po)
       end






