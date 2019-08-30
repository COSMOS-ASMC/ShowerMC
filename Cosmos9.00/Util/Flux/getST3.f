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
      real*8 ans, sum, dazim
      logical deg

      call creadParam(5)
      call cbeginRun
      call cprintPrim(ErrorOut)
      cosmax = IMAG_P(CosZenith)
      cosmin = real(CosZenith)
      azmmax = IMAG_P(Azimuth) + XaxisFromSouth
      azmmin = real(Azimuth) + XaxisFromSouth
      dazim =abs(azmmax-azmmin)*Torad
      if(cosmax .ne. cosmin) then
         write(*, *) ' this program is to compute integral of'
         write(*, *) ' primary for a fixed zenith angle '
         write(*, *) ' give the same lower and upper limit for'
         write(*, *) ' CosZenith'
         stop 9999
      elseif(  (abs(IMAG_P(Azimuth)-real(Azimuth))/360.d0 -1.d0)
     *         .gt. 1.d-3 ) then
         if(azmmax .ne. azmmin) then
            write(*,*) ' warning: fai region is not 2pi '
            write(*,*)
     *      ' integral will be performed in the given fai region'
         else
            write(*,*) ' no fai region; for a fixed fai'
            write(*,*) ' use getST2'
            stop 9999
         endif
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


      write(*,*)
     *    '          cos value=',zen1
      write(*,*)
     *    '          fai region=', real(Azimuth),' to ', 
     *               IMAG_P(Azimuth), ' deg'
      write(*,*) ' X-axis From South=', XaxisFromSouth, ' deg'

      write(*,*)
      if(rigc .eq. 0) then
         write(*,*) ' No rigidity cut is assumed'
         write(*,*)
     "   ' primary  Int(dI/dE)  sum(cumlative);'
      else
         write(*,*) ' Rigidity cut is  applied'
         write(*,*)
     "   ' primary  Int(dI/dE*RigCut)  sum(cumlative);'
      endif
      sum = 0.
      do i = 1, Prim%no_of_comps
         call inteflux(Prim%each(i), ans)
         sum = sum +  ans 
         write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum)/dazim
      enddo

        write(*,*)
     * 'If N primaries are generated in simulation,',
     * 'ST2pidcos= N/sum'
        write(*,*) 'if azimuthal angle range is not 0~2pi',
     *    ' 2pi in ST2pidcos has been adjusted'

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
         if(abs(ans/comp%inte_value-1.d0) .gt. 1.d-3 ) then
            write(*,*) ' ans=',ans, ' internal integral=',
     *           comp%inte_value
            stop 9999
         endif
      else
         call cmkptc(comp%code, comp%subcode, comp%charge, aPtcl)
         Eth =sqrt( (rigc*comp%charge)**2 + aPtcl%mass**2 )
         call cconv_prim_e2(comp, Eth, e_or_p)
         call kdwhereis(e_or_p, comp%no_of_seg+1, comp%energy, 1, j)
         if(j .le. imax) then
!              energy integral; 
            call kdexpIntFb(primdN, comp%energy(1), comp%energy(j+1),
     *           eps,  ans, error,  icon)
            ans = ans * Torad
         endif
!             add E> comp.energy(j+1)
         if(j+1 .lt. imax ) then
            call intePrim2(comp, j+1, imax, ans2)
            ans2 = ans2 * abs(azmmax - azmmin)*Torad
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


      integer icon

      real*8 eps, ans, error
      real*8 funcazim
      external funcazim
      data  eps/1.d-5/

      E = eorp
      zen = zen1
      call kdexpIntF(funcazim, azmmin, azmmax, eps, ans, error, icon)
      primdN = ans
      end
!.......................
      real*8 function funcazim(azm)
      implicit none
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zprimary.h"
      include 'Zflux.h'

      type (ptcl):: aPtcl
      real*8 rig, azm, prob, flux, zeny

!          E= E_or_P  must be converted to rigidity
      call cconv_prim_e(compx, E, aPtcl)

      rig =sqrt( aPtcl%fm%p(4)**2  - aPtcl%mass**2)/aPtcl%charge
      if(degree) then
         zeny = zen/Torad
      else
         zeny = zen
      endif
      call crigCut(azm, zeny, rig, prob)
      call cprimFlux0(compx, E, flux)  ! E = EorP
      funcazim = flux * prob 
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

