!
!        c2bdcy, c2bdcp, c2bdc0 and test programs.
!
!     include '../../KKlib/rnd.f'
!     include 'cmkptc.f'
!     include '../../KKlib/kcossn.f'
!     include 'cpxyzp.f'
!     include 'clorep.f'
!     include 'cgetRotMat4.f'
!     include 'clorez.f'
!     include 'cbst1.f'
!     include 'cbst0.f'
!     include 'cibst1.f'
!     include 'cgetlf.f'
!---------------------
!     implicit none
!:       c2bdcy, c2bdcp, c2bdc0.
!     include '../Zptcl.h'
!     include '../Zcode.h'
!     type(ptcl):: p, p1, p2
!     real*8 pabs
!
!c            test c2bdcy
!     call cmkptc(krho, 0, 0,  p)
!     p.fm.p(1)=0.
!     p.fm.p(2)=0.
!     p.fm.p(3)=10.
!     call cpxyzp(p.fm, pabs)
!     p.fm.p(4) = sqrt(p.mass**2 + pabs**2)
!
!     call cmkptc(kpion, 0,  1,  p1)
!     call cmkptc(kpion, 0, -1,  p2)
!     call c2bdcy(p, p1, p2)
!     write(*,*) p1.fm.p(1) + p2.fm.p(1)
!     write(*,*) p1.fm.p(2) + p2.fm.p(2)
!     write(*,*) p1.fm.p(3) + p2.fm.p(3)
!     call cbst1(1, p, p2, p2)
!     call cbst1(2, p, p1, p1)
!     write(*, *) p1.fm.p(3) + p2.fm.p(3)
!     write(*, *) p1.fm.p(4) + p2.fm.p(4), p.mass
!     end
!    ********************************************************
      subroutine c2bdcy(p, p1, p2)
!    ********************************************************
!        two body decay by pure phase space.  
!      Note that in case polarization exists,
!     this should not be used.
!     p: /ptcl/ Input.  particle which decays.
!               4 momentum and mass must be given
!               in a certain system K.
!    p1: /ptcl/ Input/Output. Decay product 1.
!               Mass must be given in input.
!               As output, 4 mementum is given in the
!               system K.  Other attributes are unchanged.
!    p2: /ptcl/ the same as above for 2nd ptcl.
! **** note:    for massless ptcl decay, use c2bdc0, which
!            is faster.
!
      implicit none
!----       include '../Zptcl.h'
#include  "Zptcl.h"
       type(ptcl):: p, p1, p2
       real*8 am, am1, am2, q1, e1, e2, cst, cs, sn, snt
!
        am = p%mass
        am1 = p1%mass
        am2 = p2%mass
        q1=max(  (am**2- (am1+am2)**2)*(am**2-(am1-am2)**2), 0.d0)
        q1=sqrt( q1 )/am/2
        e1=sqrt(q1**2+am1**2)
        if(am1 .eq. am2) then
           e2=e1
        else
           e2=sqrt(q1**2+am2**2)
        endif
        call rndc(cst)
        cst=cst*2-1.d0
        call kcossn(cs, sn)
        snt=sqrt(1.d0 - cst**2)
        p1%fm%p(1) = q1*snt*cs
        p1%fm%p(2) = q1*snt*sn
        p1%fm%p(4) = e1
        p1%fm%p(3) = q1*cst
        p2%fm%p(1) = -p1%fm%p(1)
        p2%fm%p(2) = -p1%fm%p(2)
        p2%fm%p(4) = e2
        p2%fm%p(3) = -p1%fm%p(3)
!         boost p1 into K
        call cibst1(1, p, p1, p1)

!          boost p2 into K
        call cibst1(2, p, p2, p2)
      end
!-----*********************************************
!:        c2bdcp test.
!     implicit none
!     include '../Zptcl.h'
!     include '../Zcode.h'
!     type(ptcl):: p, p1, p2
!     real*8 pabs
!
!     call cmkptc(kpion, 0, 0,  p)
!     p.fm.p(1)=0.
!     p.fm.p(2)=0.
!     p.fm.p(3)=10.
!     call cpxyzp(p.fm, pabs)
!     p.fm.p(4) = sqrt(p.mass**2 + pabs**2)
!
!     call cmkptc(kphoton, 0, 0,  p1)
!     call cmkptc(kphoton, 0, 0,  p2)
!     call c2bdcp(p, p1, 0.5d0, p2)
!     write(*,*) p1.fm.p(1) + p2.fm.p(1)
!     write(*,*) p1.fm.p(2) + p2.fm.p(2)
!     write(*,*) p1.fm.p(3) + p2.fm.p(3)
!       reboost to decay particle system
!     call cbst1(1, p, p2, p2)
!     call cbst1(2, p, p1, p1)
!     write(*, *) p1.fm.p(3) + p2.fm.p(3)
!     write(*, *) p1.fm.p(4) + p2.fm.p(4), p.mass
!     write(*, *) p1.fm.p(3), p2.fm.p(3), p.mass/2
!     end
!     ************************************************ 
      subroutine c2bdcp(p, p1, cst, p2)
!      2body decay when the one particle decay angle is
!     fixed.
!     ************************************************ 
!      p: /ptcl/ input.  a particle which decays. 4momentum
!           must be given in a certain system K.  
!           The mass should be given too.
!      p1:  /ptcl/.input/outut. first particle. mass must be given.
!            as input.
!            as output, 4 momentum is given at K. Others unchaed
!      cst: real*8. input.  cos(of angle to the z axis)
!             of 1st ptcl at the rest system of p
!      p2:  /ptcl/ input/output. 2nd particle. same as 1st ptcl.
!
      implicit none
!----      include '../Zptcl.h'
#include  "Zptcl.h"
      type(ptcl):: p, p1, p2
      real*8 cst
!
      real*8 am, am1, am2, q1, e1, e2, cs, sn, snt
!
      am1 = p1%mass
      am2 = p2%mass
      am = p%mass
      q1=max(  (am**2- (am1+am2)**2)*(am**2-(am1-am2)**2), 0.d0)
      q1=sqrt( q1 )/am/2
      e1=sqrt(q1**2+am1**2)
      if(am1 .eq. am2) then
           e2=e1
      else
           e2=sqrt(q1**2+am2**2)
      endif
      call kcossn(cs, sn)
      snt=sqrt(1.d0-cst**2)
      p1%fm%p(1) = q1*snt*cs
      p1%fm%p(2) = q1*snt*sn
      p1%fm%p(3) = q1*cst
      p1%fm%p(4) = e1
!
      p2%fm%p(1) = -p1%fm%p(1)
      p2%fm%p(2) = -p1%fm%p(2)
      p2%fm%p(3) = -p1%fm%p(3)
      p2%fm%p(4) = e2
!         boost p1 into K
      call cibst1(1, p, p1, p1)

!        boost p2 into K
      call cibst1(2, p, p2, p2)
      end
!-----*********************************************
!     implicit none
!:        c2bdc0 test.
!     include '../Zptcl.h'
!     include '../Zcode.h'
!     type(ptcl):: p, p1, p2
!     real*8 pabs
!
!     call cmkptc(kpion, 0, 0,  p)
!     p.fm.p(1)=0.
!     p.fm.p(2)=0.
!     p.fm.p(3)=10.
!     call cpxyzp(p.fm, pabs)
!     p.fm.p(4) = sqrt(p.mass**2 + pabs**2)
!
!     call cmkptc(kphoton, 0, 0,  p1)
!     call cmkptc(kphoton, 0, 0, p2)
!     call c2bdc0(p, p1, p2)
!     write(*,*) p1.fm.p(1) + p2.fm.p(1)
!     write(*,*) p1.fm.p(2) + p2.fm.p(2)
!     write(*,*) p1.fm.p(3) + p2.fm.p(3)
!     call cbst1(1, p, p2, p2)
!     call cbst1(2, p, p1, p1)
!     write(*, *) p1.fm.p(3) + p2.fm.p(3)
!     write(*, *) p1.fm.p(4) + p2.fm.p(4), p.mass
!     end
!     ************************************************
      subroutine c2bdc0(p, p1, p2)
!           special case for masseless 2 ptcl decay
!           such as pi0-->2g.
!     ************************************************
      implicit none
!----      include '../Zptcl.h'
#include  "Zptcl.h"
      type(ptcl):: p, p1, p2
      real*8 e1, e2, q1, cst, cs, sn, snt
!

      q1=p%mass/2
      e1=q1
      e2=e1
      call rndc(cst)
      cst=cst*2-1.d0
      call kcossn(cs, sn)
      snt=sqrt(1.d0-cst**2)
      p1%fm%p(1) = q1*snt*cs
      p1%fm%p(2) = q1*snt*sn
      p1%fm%p(4) = e1
      p1%fm%p(3) = q1*cst
!
      p2%fm%p(1) = -p1%fm%p(1)
      p2%fm%p(2) = -p1%fm%p(2)
      p2%fm%p(4) = e2
      p2%fm%p(3) = -p1%fm%p(3)
      call cibst1(1, p, p1, p1)
      call cibst1(2, p, p2, p2)
      end








