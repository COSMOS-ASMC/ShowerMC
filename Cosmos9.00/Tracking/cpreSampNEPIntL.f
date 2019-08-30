!     ******************************************************************
!     *                                                                *
!     * cpreSampNEPIntL: samples integraction length for a given
!     *              non e, gamma  particle in a given medium
!     *                                                                *
!     ******************************************************************
!     subroutine  csampNEPIntL
      subroutine  cpreSampNEPIntL(media)
      implicit none
#include  "Zglobalc.h"
#include  "Zcode.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackv.h"
#include  "Zevhnp.h"
#include  "Zmedia.h"
      type(epmedia),intent(in):: media

!     **************************************************
!
      real*8 mfp, xs, length

!            
!        m.f.p (kg/m**2) = abn /xsec(mb)
#if defined (LOOPCHK)
      integer ksave/0/, esave/-1./, ncount/0/
      integer ka
      save
      ka = TrackBefMove%p%code

      if(ksave .eq. ka .and. esave .eq.  
     *   TrackBefMove%p%fm%p(4)) then
         ncount = ncount +1
         if(ncount .gt. 10) then
            write(0,*)'  ncount=',ncount, 
     *        ' ka=',ka, ' e=',esave, ' mass=',
     *      TrackBefMove%p%mass
         endif
      else
         ksave =ka
         esave = TrackBefMove%p%fm%p(4)
         ncount =0
      endif
!///////////////////
#endif

      call cepSampNEPIntL(media, TrackBefMove%p)
      end
!!!
!!!
!!!
!!!      subroutine cknockonH
!!!      implicit none
!!!#include  "Zglobalc.h"
!!!#include  "Zcode.h"
!!!#include  "Ztrack.h"
!!!#include  "Ztrackv.h"
!!!#include  "Zelemagp.h"
!!!#include  "Zevhnv.h"
!!!      real*8 prob, path
!!!!              knock on by non e+/e- charged ptcl
!!!      call cKnockp(TrackBefMove%p, prob, path) ! path in r%l
!!!      if(prob .gt. 0.d0) then
!!!         call csetIntInf(path *X0, .false., 'knock')
!!!      else
!!!         call csetIntInf(Infty,   .false., 'knock')
!!!      endif
!!      end


      


