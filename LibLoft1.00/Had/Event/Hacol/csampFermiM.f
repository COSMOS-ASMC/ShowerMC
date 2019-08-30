      subroutine csampFermiM(t)
      implicit none
!----      include  '../../Zmass.h'
#include  "Zmass.h"
!----      include  '../../Zptcl.h'
#include  "Zptcl.h"
!
      type(fmom):: t
      real*8  ef
      real*8 u, ek, p, cs, sn, cst, snt
!           sample fermi momentum and set 4 momentum in t 
!           Fermi(surface) momentum is 200 MeV
         parameter( ef=200.e-3**2/2/masp)
!
         call  rndc(u)
         ek=ef*u
         p=sqrt(ek*(ek+masp*2))
         call kcossn(cs, sn)
         call rndc(cst)
         cst=cst*2-1.
         snt=sqrt(1.-cst**2)
         t%p(1) = p*snt*cs
         t%p(2) = p*snt*sn
         t%p(3) = p*cst
         t%p(4) = ek+masp
       end

       
