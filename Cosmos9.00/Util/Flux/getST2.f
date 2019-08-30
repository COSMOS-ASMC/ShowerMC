!          This is to compute integral of primaries for a given
!          theta, fai.
!
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
      include 'Zflux.h'

      integer  i

      real*8 cosmin, cosmax
      real*8 ans, sum
      logical deg

      call creadParam(5)
      call cbeginRun
      call cprintPrim(ErrorOut)
      cosmax = IMAG_P(CosZenith)
      cosmin = real(CosZenith)
      azmmax = IMAG_P(Azimuth) + XaxisFromSouth
      azmmin = real(Azimuth) + XaxisFromSouth

      if(cosmax .ne. cosmin .or. ( cosmax .ne. 1.d0 .and. 
     *          azmmax .ne. azmmin) ) then
         write(*,*) ' this progam is to compute integral of'
         write(*,*) ' primary for a given cos and fai'
         write(*,*) ' give the same value for upper and lower'
         write(*,*) ' limit of  CosZenith and Azimuth'
         stop
      endif
      if(cosmax .eq. 1.d0) then
         write(*,*)
     *   ' getST3 with fai=0~2pi may give you a different result',
     *   ' for cos=1'
         write(*,*) 
     *   ' both should give the same result but due to the'
         write(*,*)
     *   ' table usage, the results are different. '
      endif
      if(CutOffFile .eq. ' ') then
         rigc = 0.
         deg = .false.
      else
         rigc = 100.
         if(ZenValue .eq. 'deg') then
            deg = .true.
         else
            deg = .false.
         endif
      endif
      zen1 = cosmin
      zen2 = cosmax


      write(*,*)
      write(*,*)
     *    '          cos value=',zen1
      if(zen1 .ne. 1.d0) then
         write(*,*)
     *    '          fai value=', IMAG_P(Azimuth)
      endif
      write(*,*) '     X-axis From South=', XaxisFromSouth, ' deg'
      write(*,*)
      if(rigc .eq. 0) then
         write(*,*) ' No rigidity cut is assumed'
         write(*,*)
     *   ' primary  Int(dI/dE)  sum(cumlative);'
      else
         write(*,*) ' Rigidity cut is  applied'
         write(*,*)
     *   ' primary  Int(dI/dE*RigCut)  sum(cumlative);'
      endif
      sum = 0.
      do i = 1, Prim%no_of_comps
         call inteflux(Prim%each(i), ans)
         sum = sum +  ans 
         write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum)
      enddo

      write(*,*)     
     * 'If N primaries are generated in simulation,'
      write(*,*)' STdOmega= N/sum'

      end
!     ******************
      subroutine inteflux(comp, ans)
      implicit none
#include "Zglobalc.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
      type (component):: comp
      real*8 ans
      include 'Zflux.h'
      
      real*8 eps, error, ans2, Eth, e_or_p
      integer icon


      integer imax
      external primdN
      real*8 primdN
      type (ptcl)::aPtcl
      integer j
      data eps/1.d-4/



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
!              energy integral; 
            call kdexpIntFb(primdN, comp%energy(1), comp%energy(j+1),
     *           eps,  ans, error,  icon)
            ans = ans 
         endif
!             add E> comp.energy(j+1)
         if(j+1 .lt. imax ) then
            call intePrim2(comp, j+1, imax, ans2)
         else
            ans2 = 0.
         endif
         ans = ans + ans2
      endif
      end
    
!     ****************************      
      real*8 function primdN(eorp)
      implicit none

#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zprimary.h"
!          primary flux at E       
      real*8 eorp

      include 'Zflux.h'



      type (ptcl):: aPtcl
      real*8 rig,  prob, flux, zeny

      E = eorp
!          E= E_or_P  must be converted to rigidity
      call cconv_prim_e(compx, E, aPtcl)

      rig =sqrt( aPtcl%fm%p(4)**2  - aPtcl%mass**2)/aPtcl%charge
      if(degree) then
         zeny = zen1/Torad
      else
         zeny = zen1
      endif
      call crigCut(azmmin, zeny, rig, prob)
      call cprimFlux0(compx, E, flux)  ! E = EorP

      primdN = flux * prob 
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



