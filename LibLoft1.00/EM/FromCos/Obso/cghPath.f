!        cghPath.  mean free path of photo production of hadrons.
      subroutine  cghPath(egin, t)
!           get hadronic interaction (v.d) m.f.p for gamma
!           t is given in kg/m2
      implicit none
!

#include  "Ztrack.h"
#include  "Ztrackv.h"

      real*8 egin, t
!
      real*8  abogn, abn, xs, u
      parameter(abogn=6.02e23, abn=1.e28/abogn)  ! abn * A/ xs(mb) => kg/m2
!           gA x-section,   xs in mb
      call cgpXsec(TargetMassN, egin, xs)

      if(xs  .gt. 0.) then
         t=abn*TargetMassN /xs
         call rndc(u)
         t=-log(u)*t
      else
         t=1.0d50
      endif
      end
