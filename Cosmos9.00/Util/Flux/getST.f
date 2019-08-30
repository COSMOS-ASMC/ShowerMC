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

      real*8 cosmin, cosmax, cylr, cylh
      real*8 ans, sum1, sum2, sum3
      logical deg
      call cerrorMsg('**********IMPORTANT**********',1)
      call cerrorMsg(
     *  'You must write radius and height of a cylinder',1)
      call cerrorMsg('or 0 0 after namelist parameter',1)
      call cerrorMsg(
     *  '(1 space line may be needed before data)', 1)
      call creadParam(5)
      call cbeginRun
      call cprintPrim(ErrorOut)
      cosmax = IMAG_P(CosZenith)
      cosmin = real(CosZenith)
      azmmax = IMAG_P(Azimuth) + XaxisFromSouth
      azmmin = real(Azimuth) + XaxisFromSouth

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

      call csetCosdeg(0, .false.)
      write(*,*)
      write(*,*)
     *    '          cos region=',zen1, ' to ', zen2
      write(*,*)
     *    '          fai region=', real(Azimuth),' to ', 
     *               IMAG_P(Azimuth), ' deg'
      write(*,*) ' X-axis From South=', XaxisFromSouth, ' deg'

      write(*,*)
      if(rigc .eq. 0) then
         write(*,*) ' No rigidity cut is assumed'
         write(*,*)
     "   ' primary  Int(dI/dE)  sum(cumlative);',
     *   ' sphere'
      else
         write(*,*) ' Rigidity cut is  applied'
         write(*,*)
     "   ' primary  Int(dI/dE*RigCut)  sum1(cumlative);',
     *   ' sphere'
      endif
      sum1 = 0.
      do i = 1, Prim%no_of_comps
         call inteflux(Prim%each(i), ans)
         sum1 = sum1 +  ans 
         write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum1)
      enddo

!      ---------------------
         

      call csetCosdeg(1, deg)
      if(rigc .eq. 0) then
         write(*,*)
     *  ' primary  Int(cos x dI/dE)  sum2(cumlative)',
     *  ' sum2/sum1 flat '
      else
         write(*,*)
     *   ' primary  Int(cos x dI/dE*RigCut) sum2(cumlative)',
     *   ' sum2/sum1 flat'
      endif
      sum2 = 0.
      do i = 1, Prim%no_of_comps
         call inteflux(Prim%each(i), ans)
         sum2 = sum2 +  ans 
         write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum2), sngl(sum2/sum1)
      enddo


      call csetCosdeg(2, deg)
      if(rigc .eq. 0) then
         write(*,*)
     *  ' primary  Int((1+cos)/2 dI/dE)  sum3(cumlative);',
     *  ' sum3/sum1 hemisphere'
      else
         write(*,*)
     *  ' primary  Int((1+cos)/2 dI/dE*Rigcut) sum3(cumlative);',
     *  ' sum3/sum1  hemisphere'
      endif
      sum3 = 0.
      do i = 1, Prim%no_of_comps
         call inteflux(Prim%each(i), ans)
         sum3 = sum3 +  ans 
         write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum3), sngl(sum3/sum1)
      enddo
      write(*,*)
      read(*,*) cylr, cylh
      if(cylr .gt. 0. .and. cylh .gt. 0.) then
         ratio =2*cylr*cylh/(3.141592*cylr**2)
         call csetCosdeg(3, deg)
         if(rigc .eq. 0) then
            write(*,*)
     *      ' primary  Int((cos+ratio*sin)/sqrt(1+ratio**2)dI/dE)',
     *      '  sum3(cumlative); sum3/sum1 cylinder'
         else
            write(*,*) ' cylinder: 2rh/pir^2=ratio=',ratio
            write(*,*)
     *    ' primary  Int((cos+ratio*sin)/sqrt(1+ratio**2)dI/dE*Rc)',
     *    ' sum3(cumlative);',' sum3/sum1  cylinder'
         endif
         sum3 = 0.
         do i = 1, Prim%no_of_comps
            call inteflux(Prim%each(i), ans)
            sum3 = sum3 +  ans 
            write(*,*)' ', Prim%each(i)%symb,'  ', sngl(ans),
     *          '  ', sngl(sum3), sngl(sum3/sum1)
         enddo
      endif
      write(*,*)     
     * 'If N primaries are generated in simulation, ST= N/sum1'

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
         if(cosfactor .eq. 0) then
           if(abs(ans/comp%inte_value-1.d0) .gt. 1.d-3 ) then
              write(*,*) ' ans=',ans, ' internal integral=',
     *                    comp%inte_value
              stop 9999
           endif
           ans = ans * (zen2-zen1)* abs(azmmax - azmmin)*Torad
         elseif(cosfactor .eq. 1) then
!             take into account the horizontal  area.
            ans = comp%inte_value * (zen2**2- zen1**2)/2 *
     *                   abs(azmmax - azmmin)*Torad
         elseif(cosfactor .eq. 2) then
!              
            ans = comp%inte_value
     *           *(zen2-zen1)/2.d0* (1.d0 + (zen1 + zen2)/2)*
     *           abs(azmmax - azmmin)*Torad
         elseif(cosfactor .eq. 3) then
!                  cos + rat*sin
            ans = comp%inte_value *
     *       (  (zen2-zen1)/2.0d0 + ratio* (
     *           acos(zen1)-acos(zen2) +
     *           0.5*( zen2*sqrt(1.-zen2**2) -
     *                 zen1*sqrt(1.-zen1**2) )) )
     *           /sqrt(1.+ratio**2)
     *          * abs(azmmax - azmmin)*Torad
         else 
            write(*,*) ' error of cosin'
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
            if(cosfactor .eq. 0 ) then
               ans2 = ans2 * (zen2- zen1)*abs(azmmax - azmmin)*Torad
            elseif(cosfactor .eq. 1) then
               ans2 = ans2 * (zen2**2- zen1**2)/2 *
     *                  abs(azmmax - azmmin)*Torad
            elseif(cosfactor .eq. 2) then
!              fai is from 0 to 2pi; 
               ans2 =
     *          ans2 *(zen2-zen1)/2.d0* (1.d0 + (zen1 + zen2)/2.d0)
            else
                ans2 = ans2 *
     *          (  (zen2-zen1)/2.0d0 + ratio* (
     *           acos(zen1)-acos(zen2) +
     *           0.5*( zen2*sqrt(1.-zen2**2) -
     *                 zen1*sqrt(1.-zen1**2) )) )
     *           /sqrt(1.+ratio**2)
     *           *  abs(azmmax - azmmin)*Torad
            endif
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

      real*8 funczen, ans
      external funczen

      E = eorp
!         integral on  zenith 
      call k16pGaussLeg(funczen, zen1, zen2, 16, ans)
      primdN = ans
      end
!      ********************************
      real*8 function funczen(zenx)
      implicit none
#include "Zptcl.h"
#include "Zprimary.h"

      real*8 zenx

      include 'Zflux.h'

      integer icon
      real*8 eps, ans, error
      real*8 funcazim
      external funcazim
      data  eps/1.d-5/

      zen = zenx
      call kdexpIntF(funcazim, azmmin, azmmax, eps, ans, error, icon)
      funczen = ans
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
!/////////////
!      prob = 1.
!/////////////
      call cprimFlux0(compx, E, flux)  ! E = EorP
      if(cosfactor .eq. 0 ) then
         funcazim = flux * prob 
      elseif(cosfactor .eq. 1) then
         funcazim = flux * prob * zen
      elseif(cosfactor .eq. 2) then
         funcazim = flux * prob *(1.0d0+ zen)/2.d0
      else
         funcazim = flux * prob *
     *   (zen + ratio*sqrt(1.-zen**2)) / sqrt(1.+ratio**2)
      endif
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
!     ***************************
      subroutine csetCosdeg(cosin,  degin)
      implicit none
      logical cosin, degin
#include "Zptcl.h"
#include "Zprimary.h"

      include 'Zflux.h'
      
      cosfactor = cosin
      degree = degin
      end


