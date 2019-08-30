!     ****************************************
!     *                                                          *
!     * cmulScat: multiple Coulomb scattering
!     *                                                          *
!     ****************************************
!
!
      subroutine cmulScat(media, theta)
!      use modXsecMedia
      use modMCScontrol
      use modSetIntInf
      implicit none
!       Using  cTrack and Move, compute scattering angle 
!
#include "Zglobalc.h"
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
#include  "Zmedia.h"      


      type(epmedia),intent(in):: media
      real*8 theta         ! output. sampled angle in radian.
!
!            
!
      integer cond
      real(8):: sgr
      integer ic

!             sample theta
      if(IntInfArray(ProcessNo)%length .gt. 0.) then
         if(ActiveMCS == "El_con") then
            !  total length in kg/cm2; sample theta                   
            sgr =IntInfArray(ProcessNo)%thickness
            call cElsepaCondense(sgr, theta)
         elseif(abs(Moliere) >= 2)  then   ! abs is historical reason
            call cSampMol(media, theta, cond)  ! rigorous Moliere
            if(cond /= 0 ) then
               call cmulScat2(theta)   ! Gauss
            endif
         elseif(abs(Moliere) == 1)  then
            call cmulScat1(theta)      ! simplified Moliere
         elseif( Moliere == 0 ) then
            call cmulScat2(theta)      ! Gauss
         endif
      else
         theta = 0.
      endif
      end
!     **************************
      subroutine cmulScat2(media, theta)
      use modSetIntInf
      implicit none
!       Gaussian scattring
!
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
#include  "Zmedia.h"
      
      type(epmedia),intent(in):: media

      real*8 theta ! output. sampled spatial angle in radain.

      real*8 tetarms, g1, g2, u, beta2, dt
      integer nc
      
      real*8 hpi 
      parameter(hpi = pi/2.)

      g1 = TrackBefMove%p%fm%p(4)/MovedTrack%p%mass
      g2 = MovedTrack%p%fm%p(4)/MovedTrack%p%mass
      beta2 = 1.d0 - 1.d0/g1/g2
      if(beta2 .le. 0.) then
         tetarms = 0.
      else
         dt = IntInfArray(ProcessNo)%thickness/ media%X0  ! dt  r%l
         if(dt .gt. 1.d-3) then
            tetarms = Es/TrackBefMove%p%fm%p(4)*
     *       abs(MovedTrack%p%charge) *
     *       sqrt(dt/beta2)*(1.0 + 0.038*log(dt))
         else
            tetarms = Es/TrackBefMove%p%fm%p(4)*
     *      abs(MovedTrack%p%charge) * sqrt(dt/beta2)
         endif
      endif
      theta = pi
      nc = 0
      do while(theta .gt. hpi)
         if(nc .gt. 10) then
!              tetarms seems too large
            theta = u**0.1 * hpi  ! give some value 
         else
            call rndc(u)
            theta = sqrt(-log(u))* tetarms
            nc = nc +1
         endif
      enddo
      end
      subroutine cmulScat1(theta)
      use modSetIntInf
      implicit none
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
!           by Moliere thoery 

      real*8 theta
      
      real*8 g1, g2, leng
      integer cond

      real*8 t, tmp, avx, avy, disp, cs, sn, e1, e2, d1, d2
      real*8 rho, cvh2den
      integer chg
 
      e1 = TrackBefMove%p%fm%p(4)
      g1 = e1/TrackBefMove%p%mass
      e2 = MovedTrack%p%fm%p(4)
      g2 = e2/MovedTrack%p%mass
      rho = cvh2den(
     *      (TrackBefMove%pos%height+MovedTrack%pos%height)/2.d0
     *      )

      leng = IntInfArray(ProcessNo)%length  ! in m

      chg =  MovedTrack%p%charge
      call cmoliere(rho,  chg, MovedTrack%p%mass, g1, g2,
     *     leng, theta, cond)
      if(cond .ne. 0) then
!           Moliere theory cannot be applied
         call cmulScat2(theta)
      endif
      end

      subroutine cElsepaCondense(sgr, theta)
!     can be eplaced by eElsepaCondense(sgr, theta)
!         but sgr is not kg/m2 (g/cm2)      
      use modMCScontrol
      use modsoftMCS
      use modcMCS
      implicit none
#include "ZmediaLoft.h"
      
      real(8),intent(in):: sgr  ! path length in kg/m2
      real(8),intent(out):: theta ! sampled angle
      
      real(8)::cm2topgrm,s0,s1, s2, ls1, ls2 
      real(8):: mu, cosa, sina
!!!      call ciniFsoft(KEeV, sgr)   is for integration
!!!    from 0 to muc.  We need here from 0 to 1. so
!!!   do below
!                         1 is to get s in cm2
      call cTPXS(TPXSnow, 1,  KEeV, s0, s1, s2)
!                              not 0.1 (/(kg/m2))-->(/g/cm2)
      cm2topgrm = Media(MediaNo)%mbtoPkgrm/ 1d-27 
      ls1 =1.0/(s1*cm2topgrm)
      ls2 =1.0/(s2*cm2topgrm)

      avemu = (1.0d0 - exp(-sgr/ls1))/2
      avemu2 = avemu - (1- exp(-sgr/ls2) )/6.d0
      if( avemu > 0.49999999d0 ) then
          ! mu is uniform (asoft =1 bsoft =1 is also ok)
         bsoft = 0.5d0
         asoft = 0.5d0
      else
         bsoft =( 2*avemu - 3*avemu2 )/(1-2*avemu)
         asoft = (1-2*avemu) + bsoft
      endif
      call csampSoftMCSang(mu, cosa)
      theta = acos(cosa) 
      end
