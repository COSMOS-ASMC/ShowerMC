!           testing cgetcm; getting cm frame of two particles
!     include 'cgetlf.f'
!     include 'clorez.f'
!     include 'cgeqm.f'
! --------------------------------------
!
!     implicit none
!     include '../Zptcl.h'
!
!     type(fmom):: gc, p1, p2, gci
!     type(ptcl):: pp1, pp2, t1, t2
!     integer icon
!     p1.x=0.
!     p1.y=0.
!     p1.z=1.e3
!     p2.x=0.
!     p2.y=0.
!     p2.z=0.
!
!     pp1.mass=.938
!     pp2.mass=.938
!     p1.e=sqrt(p1.x**2 + p1.y**2+p1.z**2 + pp1.mass**2)
!     p2.e=sqrt(p2.x**2 + p2.y**2+p2.z**2 + pp2.mass**2)
!     pp1.fm = p1
!     pp2.fm = p2
!     call cgetcm(pp1, pp2, gc, icon)
!     write(*, *) gc.p
!     t1 = pp1
!     t2 = pp2
!     gci.x = - gc.x
!     gci.y = - gc.y
!     gci.z = - gc.z
!     gci.t = gc.t
!     call clorez(gci, pp1,  t1)
!     call clorez(gci, pp2,  t2)
!     write(*, *)  t1.fm.p, t2.fm.p
!     call cgetcm(t1, t2, gc, icon)
!     write(*, *) gc.p
!     end
      subroutine cgetcm(p1, p2,  gc, icon)
!         get Lorentz factor of the Center of Momentum
!         system of arbitrary two particles.
!      p1: type ptcl. Input.  particle 1
!      p2: type ptcl. Input.  particle 2
!      gc: type fmom. Output. (g*beta, g) of the
!          cm frame. 
!    icon: integer. output.  0 -> ok.  1 -> no cms definable.
!              if icon = 1, gc is undef.
         implicit none
!----         include '../Zptcl.h'
#include  "Zptcl.h"
         type(ptcl):: p1, p2
         type(fmom):: gc
!
         type(fmom):: q

         integer icon
!
         call cgeqm(p1, p2, q, icon)         
         if(icon .eq. 0) then
!            get Lorentz factor of q
             call cgetlf(q, gc)
         endif    
        end         

