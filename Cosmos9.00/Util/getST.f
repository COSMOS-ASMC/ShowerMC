#include  "BlockData/cblkGene.h"
      implicit none
#include "Zglobalc.h"
#include "Zmanagerp.h"
#include "ZrigCut.h"
#include "Zptcl.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zprimary.h"
#include "Zprimaryc.h"
#include "Zprimaryv.h"
#include "Zincidentp.h"
#ifdef NEXT486
#define IMAG_P dimag
#else
#define IMAG_P imag
#endif


      integer  i

      real*8 cosmin, cosmax
      real*8 ans, sum
      
      real*8  azmmin, azmmax, zen1, zen2, rigc
      type (component)::compx
      common/ZgetST/compx, zen1, zen2, azmmin, azmmax, rigc

      call creadParam(5)
      call cbeginRun
      call cprintPrim(ErrorOut)
      cosmax = IMAG_P(CosZenith)
      cosmin = real(CosZenith)
      azmmax = IMAG_P(Azimuth) + XaxisFromSouth
      azmmin = real(Azimuth) + XaxisFromSouth
      if(CutOffFile .eq. ' ') then
         rigc = 0.
         call csetCosdeg(.true., .false.)
      else
         rigc = 100.
         if(ZenValue .eq. 'deg') then
            call csetCosdeg(.true., .true.)
         else
            call csetCosdeg(.true., .false.)
         endif
      endif
      zen1 = cosmin
      zen2 = cosmax

      write(*,*) ' primary  Int(cos x dI/dE)  sum(cumlative)'
      sum = 0.
      do i = 1, Prim%no_of_comps
         call inteflux(Prim%each(i), ans)
         sum = sum +
     *     ans * (cosmax-cosmin)*abs((azmmax-azmmin)*Torad)
         write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum)
      enddo
      write(*,*)
     *  'If N primaries are generated in simulation, ST= N/sum'
      end
!     ******************
      subroutine inteflux(comp, ans)
      implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
      type (component):: comp
      real*8 ans
      
      real*8  azmmin, azmmax, zen1, zen2, rigc
      type (component)::compx
      common/ZgetST/compx, zen1, zen2, azmmin, azmmax, rigc

      
      real*8 eps, error, ans2, Eth, e_or_p
      integer icon


      integer imax
      external primdN
      real*8 primdN
      type (ptcl)::aPtcl
      integer j
      data eps/1.d-4/

!          take into account the horizontal  area.

      compx = comp
      imax = comp%no_of_seg
      if(rigc .eq. 0.) then
!         no rigidy  cut. integrate segmented power functions
         call intePrim2(comp, 1, imax, ans)
      else
         call cmkptc(comp%code, comp%subcode, comp%charge, aPtcl)
         Eth =sqrt( (rigc*comp%charge)**2 + aPtcl%mass**2 )
         call cconv_prim_e2(comp, Eth, e_or_p)
         call kdwhereis(e_or_p, comp%no_of_seg+1, comp%energy, 1, j)
         if(j .le. imax) then
            call kdexpIntFb(primdN, comp%energy(1), comp%energy(j+1),
     *           eps,  ans, error,  icon)
!             add E> comp.energy(j+1)
         endif
         if(j+1 .lt. imax ) then
            call intePrim2(comp, j+1, imax, ans2)
         else
            ans2 = 0.
         endif
         ans = ans + ans2
      endif
      end
    
!     ****************************      
      real*8 function primdN(E)
      implicit none

#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zprimary.h"
!          primary flux at E       
      real*8 E

      real*8  azmmin, azmmax, zen1, zen2, rigc
      type (component)::compx
      common/ZgetST/compx, zen1, zen2, azmmin, azmmax, rigc

      real*8 flux
      call cprimFlux(compx, E, rigc, zen1, zen2, azmmin,
     *     azmmax, flux)
      primdN = flux
      end
      subroutine intePrim2(comp, i1, i2, ans)
      implicit none
#include "Zptcl.h"
#include "Zprimary.h"
!          primary flux at E       
      type (component)::comp
      integer i1, i2
      real*8 ans

      real*8 sum, beta

      integer i
      sum = 0.

      do i = i1, i2
         beta=comp%beta(i)
         if(beta .ne. 1.) then
            sum = sum +
     *           comp%flux(i)*comp%energy(i)
     *           / (beta - 1.)
     *           *( 1. -
     *           (comp%energy(i+1)/comp%energy(i))**(1.-beta))
         else
            sum =
     *           sum +  comp%flux(i)* comp%energy(i)
     *           * log(comp%energy(i+1)/comp%energy(i))
         endif
      enddo
      ans = sum
      end
      subroutine chookTrace
      end
      subroutine chookCeren
      end
      subroutine chookCerenS
      end
      subroutine chookCerenE
      end
      subroutine chookBgRun
      end
