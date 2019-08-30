      subroutine epSampMol(mediax, Mol, teta, cond)
!      use modMultipleScat       not needed
      implicit none
!        Moliere theory of multiple scattering angle.
!        with improvement by Bethe.  rejection factor sqrt(sin teta/teta)
!        is included 
!        
!        The other is the next factor which need not be applied ?
!         exp(Bxc2/16.) rejection may be tried. We use Moliere when Bxc2 <1
!         So, if this rejection is to be tried, 
!         since for Bxc2=1, exp(..)=1.064 so the rejection for
!         smaller Bxc2 must be done by
!         u <  exp(Bxc2/16)/1.064. 
!      Define next for this correction;
!cccc # define REJBYEXPBXC2
!               For Pb, 1GeV e, there is no difference with this correction
!         of without this correction.
!     
#include "Zglobalc.h"
#include "Zmass.h"
#include "Zmedia.h"      
! #include "ZepTrackv.h"

      type(epmedia):: mediax    ! input. media 
      integer,intent(in):: Mol  ! same as Moliere.
                           
      real(8),intent(out)::teta !  sampled spatial angle in radian.
      integer,intent(out)::cond !  0 ok. non-0. Moliere theory 
       ! not applicable. Also, if Mol(iere)=4 and exp(-x) sampling
       ! is accepted, no actual sampling is done here and cond =2
       ! is set. (Later, Gaussin should be tried).
      integer::icon
      real(8):: u
      logical ok

      real*8  B, xc2, Bxc2
      common /Zmoliere/ B, xc2, Bxc2

#ifdef REJBYEXPBXC2
      real(8):: expBxc2
#endif      
      ok = .false.
      cond = 0
      do while (.not. ok)
         call epSampMol0(mediax, Mol, teta, icon)
         if( icon == 0 ) then !  Moliere
            call rndc(u)
!                
#ifdef REJBYEXPBXC2
            expBxc2 = exp(Bxc2/16.0d0)/1.064d0  !   needed???
#endif
!
            if(teta < 0.5d0 ) then 
!                       sin(teta) /teta = (teta - teta^3/6)/teta
!                      = 1- teta^2/6.  max error  0.054 %
               ok = u**2 < (1.0d0- teta**2/6.0d0)
#ifdef REJBYEXPBXC2
!                    needed ???
               if(ok) then 
                  call rndc(u)
                  ok =  u < expBxc2
               endif
#endif
            elseif( teta < 3.1415 ) then
               ok = u**2 < sin(teta)/teta
#ifdef REJBYEXPBXC2
               if(ok) then 
                  call rndc(u)
                  ok =  u < expBxc2
               endif
#endif
            else
               ok = .false.
            endif
         elseif( icon == 1) then  !  GS
!             rejeciton not needed; GS is used
            ok =.true.
         else  ! should use Gaussian
            ok =.true.
            cond = 1
         endif
      enddo
      end
      
      
      subroutine epSampMol0(mediax, Mol, teta, cond)
      use  SampMolReducedA
      use  modMultipleScat
      implicit none
!        Moliere theory of multiple scattering angle.
!        with improvement by Bethe.
#include "Zglobalc.h"
#include "Zmass.h"
#include "Zmedia.h"
! #include "ZepTrackv.h"

      type(epmedia):: mediax    ! input. media 
      integer,intent(in):: Mol  ! same as Moliere.
      real(8),intent(out)::teta !  sampled spatial angle in radian.
      integer,intent(out)::cond !  0 ok.  Moliere used
                                !  1 ok.  G.S is used.
                                !  2 n.g  Gauss should be used.

      real*8  B, xc2, Bxc2
      common /Zmoliere/ B, xc2, Bxc2


      
!     *********************
      real*8  expb
      real(8):: sb,  tetasq

      call epPreScatCalc  ! get g*beta^2, beta^2, (Me/M)^2 in Zcmdeia
!          get Xc^2
      call epXaic2(mediax,   xc2)
!           we don't get xa2 but directly get exp(b)
      call epMolExpb(mediax,   expb)
      sb =log(expb)  ! 
!///////////
!      write(*,'(a, 1p,6g14.5)')
!     *  ' sb ', sb, beta2, Move.dt, Move.dl, cTrack.p.fm.p(4),
!     *          Move.track.p.fm.p(4)
!//////////////
!          B -log(B) = b  
!///      if(sb .lt. 3.395d0) then
!!!       if(sb .lt. 1.5d0) then
      if(  sb < sbMin ) then    !  defined in module  SampMolReducedA
!         Moliere theory cannot be appliled; use Gaussian later
!         (almost no scattering)
         cond = 2
      else
         call epMolBlogB(sb, B)  !  B
         Bxc2 = B*xc2
         if(Bxc2 > 1.0d0) then
!             normally this happens for low E, which stops soon
!             B  so not serious.
!               use GS method ? --> cond=1
            cond =2   ! use Gauss at present
         else
            call epSampMolReducedA(B, Mol, tetasq, cond)
            if( cond == 0 ) then   
               teta = sqrt(tetasq *Bxc2)
            endif
!///////////
!      write(*,'(a, i2, 1p,9g14.5)')
!     *  'B ', cond, sb,beta2,B, Move.dt, Move.dl, cTrack.p.fm.p(4), 
!     *        Move.track.p.fm.p(4),  Bxc2, teta
!//////////////
         endif
      endif
      end 
!
      subroutine epXaic2(mediax, xc2)
!       compute xc2 of Moliere's theory of scattring.
!      epPreScatCalc must have been called beforehand for 
!      each scattering.
      use modMultipleScat
      implicit none
#include "Zglobalc.h"
#include "Zmedia.h"
#include "Zmass.h"      
! #include "ZepTrackv.h"
#include "Zcode.h"

      
      type(epmedia):: mediax    ! input. media
      real*8 xc2   ! output. Xc^2 in radian^2





!
!      integeral 0 to leng of  1/(beta**4 E**2)
!      is approximated as  leng/(beta1**2 gamma1 beta2**2 gamma2)
!      Note: bata**2 *gamma = gamma - 1/gamma
!


!     if( cTrack%p%code == kelec) then
      if( code == kelec ) then      
         xc2 = mediax%MoliereForXc2
      else
!     xc2 = 0.6011d0*cTrack%p%charge**2 * massratio2 *
         xc2 = 0.6011d0 * charge**2 * massratio2 *
     *     mediax%Z2/mediax%A 
      endif
!      xc2 = xc2*Move%dx/gbeta2  ! dx is in g/cm2. no need to consider rhoc.
! real length is dx/(rho*rhoc)
      xc2 = xc2 * dxgpcm2 /gbeta2  ! dx is in g/cm2. no need to consider rhoc.      
      end

      subroutine epMolExpb(mediax, expb)
!            compute exp(b) of Moliere's scattering theory.
!            epPreScatCalc must have been called beforehand for 
!     each scattering.
      use modMultipleScat
      implicit none
#include "Zglobalc.h"
#include "Zmedia.h"
!#include "ZepTrackv.h"

      
      type(epmedia):: mediax    ! input. media
      real*8 expb   ! output. exp(b) of Moliere's theory.

      
      real(8):: sum, alfai2, temp
      integer::i 

      if( beta2 > 0.99 ) then  ! beta > 0.995 --> regards 1
         sum = mediax%MoliereExpb
      else
         sum =0.
!     temp = (cTrack%p%charge/137.0)**2/beta2
         temp = (charge/137.0)**2/beta2
         do i = 1,  mediax%noOfElem

            alfai2 = mediax%elem(i)%Z**2 *temp

            sum = sum +
     *       mediax%No(i) *
     *       mediax%elem(i)%Z**0.3333*(mediax%elem(i)%Z + 1.d0)
     *       / (1.d0 + 3.327d0*alfai2)
         enddo
         sum = 6702.d0 *sum/mediax%A
      endif
!     expb = sum *Move%dx * cTrack%p%charge**2 /beta2
      expb = sum * dxgpcm2 * charge**2 /beta2      
      end

