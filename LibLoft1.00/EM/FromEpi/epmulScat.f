!     ****************************************
!     *                                                              *
!     * epmulScat: multiple Coulomb scattering
!     *                                                              *
!     ****************************************
!
!
      module modMultipleScat
      implicit none
      real*8  gbeta2, beta2, massratio2

      
      real(8):: Ebef
      real(8):: Eaft
      real(8):: dlcm
      real(8):: dxgpcm2
      real(8):: dtrl
      real(8):: mass
      integer:: code
      integer:: charge
      
      logical:: ready      
      end  module modMultipleScat
      
      subroutine epmulScat(Pbef, Paft, dl, dx, mediax,  theta)
      use modMCScontrol
      use modMultipleScat
      use modEMcontrol
      implicit none
!       Using  cTrack and Move, compute scattering angle 
!
#include "Zglobalc.h"
#include "Zptcl.h"
#include "Zmedia.h"
! #include "ZepTrackp.h"
!  #include "Zelemagp.h"      
      

      type(ptcl),intent(in):: Pbef, Paft ! ptcl before and after where
!                                 scattering occurs.
      real(8),intent(in):: dl !  traveling distance between before and after
      real(8),intent(in):: dx   !  dl in cm ; dx in g/cm2.  It is
!                              dl*mediax%rho*mediax%rhoc but
!                              for Cosmos not so (rho changes)
      type(epmedia),intent(in):: mediax
      real(8),intent(out):: theta ! scattering angle in rad.
!----------------------
      integer cond
      integer ic
!----------------

      
      Ebef = Pbef%fm%p(4)
      Eaft = Paft%fm%p(4)
      dlcm = dl
      dxgpcm2 = dx
      dtrl    = dx / mediax%X0g !  radiation length  ok for Cosmos too.
      code = Pbef%code
      charge = Pbef%charge
      mass = Pbef%mass
!----------------
!            ! gbeta etc not yet obtained
      ready = .false.
!!!!!!!!!!
!      write(0,*) ' mulScat: Moliere=', Moliere
!      write(0,*) ' ActiveMCS=', ActiveMCS
!!!!!!!!!!      
!
!             sample theta
      if(dlcm .gt. 0.) then
!     
         if(ActiveMCS == 'El_con') then
            call epElsepaCondense(dxgpcm2, theta)
         elseif(Moliere == 1  .or. Moliere == 3 ) then
!          Cosmos unit must be used.
!             ic = cTrack%p%charge  !  *2 --> *4 conversion;  not used
            call epmoliere(mediax, Moliere, theta, cond)
            if(cond .ne. 0) then
!               call epang2(Move.dt, 
!     *              cTrack.p.fm.p(4)*Move.Track.p.fm.p(4),
!     *              theta)
               call epmulGauss(theta)
            endif
         elseif( Moliere >=  2) then
            call epSampMol(mediax, Moliere, theta, cond)
            if(cond .ne. 0) then
               call epmulGauss(theta)
            endif
         else
            call epmulGauss(theta)
         endif
      else
         theta = 0.
      endif
      end
!!          not used now      
!!      subroutine epang2(t, e2, teta)
!!      implicit none
!!!              simple  gaussian  approx. obsolute but
!!!         usable.
!!#include "ZepTrackp.h"
!!      real*8 t  ! input. path length in r.l
!!      real*8 e2 ! input  E1 x E2  at path head and end
!!      real*8 teta ! output.  sampled angle in radian
!!
!!      real*8 u
!!      call rndc(u)
!!           
!!      teta=Escat*sqrt(max(-log(u)*t/e2, 0.d0))
!!      end
!!!     **************************
      subroutine epmulGauss(teta)
      use modMultipleScat
      use modEMcontrol
      implicit none
!       Gaussian aprox.
!
#include "Zglobalc.h"
! #include "ZepTrackp.h"
!  #include "ZepTrackv.h"

      real*8 teta ! output. sampled spatial angle in radain.

      real*8 tetarms, u
      integer nc
      real*8 hpi 
      parameter(hpi = pi/2.)

      call epPreScatCalc
!     g1 = cTrack%p%fm%p(4)/cTrack%p%mass
!      g1 = Ebef/mass
!      g2 = Eaft/mass
!      beta2 = 1.d0 - 1.d0/g1/g2   ! < beta^2>
      if(beta2 .le. 0.) then
         tetarms = 0.
      else
!     if(Move%dt .gt. 1.d-3) then
!!!!!!!!!!
!         write(0,*) ' gauss: dtrl=',dtrl, ' Ebef, Eaft=',
!     *   Ebef, Eaft
!!!!!!!!!!         
         if(dtrl .gt. 1.d-3) then
!                Es/(pbeta) = Es/(mgbeta beta)=Es/mg/beta2
!                Es/E/beta2
!     tetarms = Escat/cTrack%p%fm%p(4)*abs(cTrack%p%charge) *
!     *        sqrt(Move%dt)/beta2*(1.0 + 0.038*log(Move%dt))            
            tetarms = Escat/Ebef*abs(charge) *
     *        sqrt(dtrl)/beta2*(1.0 + 0.038*log(dtrl))  
!     *        sqrt(Move.dt)/beta2
         else
            tetarms = Escat/Ebef*
     *          abs(charge) * sqrt(dtrl)/beta2
         endif
!!!!!!!!!!!
!         write(0,*) ' tetarms =', tetarms
!!!!!!!!!1111         
      endif
      teta = pi
      nc = 0
      do while(teta .gt. hpi)
         if(nc .gt. 10) then
!              tetarms seems too large
            teta = u**0.1 * hpi  ! give some value 
         else
            call rndc(u)
            teta = sqrt(-log(u))* tetarms
            nc = nc +1
         endif
      enddo
      end
!     ***************************************
      subroutine epmoliere(mediax, Mol, teta, cond)
      use modMultipleScat
      use SampMolReducedA
      implicit none
!        Moliere theory of multiple scattering angle.
!        This is a modified version from the one used in Cosmos
!     change: interface.    teta > pi/2 dose not appear.
#include "Zglobalc.h"
#include "Zmedia.h"      
#include "Zmass.h"
!!!!!!!!!!!!!!
#include "Zptcl.h"
!!!!!!!!!!!!      
!  #include "ZepTrackv.h"

      type(epmedia):: mediax    ! input. media 
      integer,intent(in):: Mol ! same as Moliere
      real*8 teta ! output. sampled spatial angle in radin.
      integer cond ! output. 0 ok. non-0. Moliere theory 
        ! not applicable or, if Mol(iere)=3, and exp(-x)
        ! sampling is selected, no sampling is done here
        ! so that later, Gaussian sampling should be tried

      real*8 hpi2
      parameter (hpi2 = (pi/2)**2 )

      real*8  b, xc2
      common /Zmoliere/ b, xc2


      
!     *********************
      real*8 xa2, bp,  u
      real*8 a0, a1, a2,  sum, ra2, ra2inv
      real*8 g1  !  gamma factor at the path head
      real*8 g2  !  gamma factor at the path end
      integer icon
      real*8 rejf1, rejf21, rejf22
      real*8 x
!       
!      rejection function for redueced angle < 1.8
!          x is square of reduced angle better than 0.2 %
       rejf1(x)= ((0.1217176d-01*x + 0.3054916d-01)*x -0.2524543d0)*x
     *          + 0.9990352d0          
!
!              at 0 < x = 1/angle^2 < 0.15
       rejf21(x) =(( -162.1568*x + 44.48334)*x + 0.3907116)*x 
     *      + 0.4399338              

!             at  x = 1/angle^2 > 0.15
       rejf22(x) = (( 71.23819*x - 49.61906)*x + 10.77501)*x+ 0.2001145      
!

! ------------------------------------
!     g1 = cTrack%p%fm%p(4)/cTrack%p%mass
       call epPreScatCalc
!       g1 = Ebef/mass
!       g2 = Eaft/mass
!       gbeta2=(g1 - 1.d0/g1) * (g2 - 1.d0/g2) ! < (g*beta^2)^2 >
!       beta2 = 1.d0 - 1.d0/g1/g2 !        < beta^2>
!       massratio2= (masele/mass)**2 ! (me/m)^2

!    ............................

!          get Xc^2
       call epkaic2(mediax,   xc2)
!          get Xa^2
      call epkaia2(mediax,   xa2)

!          b -log(b) = b'
      bp = log(xc2/xa2/1.167)


!      if(bp .lt. 3.395) then
      if(bp .lt. sbMin) then    ! v9.201
!         Moliere theory cannot be appliled; use Gaussian later
!         (almost no scattering)
         cond = 1
      else
         cond = 0
         call epblogb(bp, b, icon)
         a0 = max(1.d0 - 5/b, 0.d0)  ! use single scattering term if b<=5. 
         icon = 1               ! make 0 if no rejection
!           the sampling function decomposition is explained in Test/....tex
!     |  Cosmos |
         do while (icon .ne. 0)
            a1 = 5.21062/b
            a2 = 0.7128/b
            sum = a0 + a1 + a2
            call rndc(u)
            if(a0/sum  .gt. u) then
!              sample reduced angle from exp(-x) dx where x = reduced
!              angle^2.
               if( Mol == 3 ) then  ! v9.170
                  cond =  2
                  return
               endif
               call rndc(u)
               ra2 = -log(u)
               icon = 0
            elseif( (a0+a1)/sum .gt. u) then
!                 sample reduced angle from exp(-x) dx (same as above but
!                in the region of ra < 1.8
               call rndc(u)
               ra2 = -log(1.-u/1.04076)

!                 rejection function
               call rndc(u)
               if(u .lt. rejf1(ra2)) then
!                 !!
                  if(Mol == 3 ) then
                     cond=3
                     return
                  endif
                  icon = 0
               endif
            else
!                 sample reduced angle from 2xc2 x^-4dx
               call rndc(u)
               ra2 = 3.24/u
!                   rejection function
               call rndc(u)
               ra2inv = 1./ra2
               if(ra2inv .lt. 0.15) then
                  if(u .lt. rejf21(ra2inv)) then
                     icon = 0
                  endif
               elseif(u .lt. rejf22(ra2inv)) then
                  icon = 0
               endif
            endif
            teta = ra2 * xc2 * b !  actually  theta ^ 2
            if(teta .ge. hpi2 ) then
               icon = 1
            else
               icon = 0
               teta = sqrt(teta)               
            endif
         enddo
      endif
      end 
!     *******************
      subroutine  epqmoliere(bb,  xxc2)
!     ******************
      implicit none
      real*8 bb, xxc2
      real*8  b, xc2
      common /Zmoliere/ b, xc2


!            inquire the latest consts for Moliere function.
!      You can make reduced angle from theta by
!          reduced angle^2 = theta^2/b/xc2
      bb = b
      xxc2 = xc2
      end

!     ***************************
      subroutine epblogb(c, b, cond)
      implicit none
!        solve  B - log(B) = c
!
      real*8 c ! input.   c>=1.
      real*8 b ! output.  solved b >=1. (b <1 is discarded)
      integer cond ! output. 0 if ok.
                   !         1 if c < 1.
      if(c .lt .1) then
         cond = 1
      else
!         b = 0.7 + 1.32 *c -  0.01451* c*c
          b =(((-0.3733404E-04*c + 0.1902303E-02)*c -0.3841290E-01 )*c
     *         + 1.431932)*c +     0.5200487 
          cond = 0
      endif
      end
!     ***************************
      subroutine epkaia2(mediax,   xa2)
      use modMultipleScat
      implicit none
#include "Zglobalc.h"
!     #include "ZepTrackv.h"
#include "Zmedia.h"      

!        compute Xa^2; assume the Xa^2 is weakly
!       dependent on Z, we use average Z=zave for
!       calculation.
!
      type(epmedia):: mediax
 
      real*8 xa2   !  output. Xa^2.

      real*8   alpha, const, large
      parameter (alpha = 1./137., const = (1.13*alpha)**2 )
      parameter (large = (pi/2.)**2)
      real*8 sum1, sum2, temp, xai2
      integer i
!/////////
!      real*8 xa2temp
!      write(*,*) ' gbeta2 ', gbeta2, massratio2
!////////////

      if(gbeta2 .le. 0.) then
         xa2 = large
!/////////
!        xa2temp = 0.
!////////////
      else
!//////////
!        This is normally very good.  For heavy compound media 
!        like PWO/BGO, 10 % diff. (larger than the method
!        below) 
!        but we don't see diff in the resultant cascade spread 
!         xa2temp  = const * mediax.Zeff**0.66666 * massratio2 *
!     *   (1.13 * beta2 + 3.76*(alpha*cTrack.p.charge*mediax.Zeff)**2)
!     *    /gbeta2
!///////////
         sum1 =0.
         sum2 = 0.

         if( mediax%noOfElem .eq. 1 )  then
            xa2 = const*mediax%Z2_3rd*massratio2 *
     *           (1.13*beta2 +
     *          3.76*(alpha*charge*mediax%Z)**2)/gbeta2
!     *          3.76*(alpha*cTrack%p%charge*mediax%Z)**2)/gbeta2
         else
            do i = 1, mediax%noOfElem
               xai2 = const*mediax%elem(i)%Z**0.666*massratio2 *
     *          (1.13*beta2+ 
     *        3.76*(alpha * charge * mediax%elem(i)%Z)**2)/ gbeta2
!                  even heavy Z(Z+1) is recommended as well as e+/e-
               temp = mediax%No(i)*
     *         (mediax%elem(i)%Z + 1.)*mediax%elem(i)%Z 
               sum1 = sum1+ temp
               sum2 = sum2 + temp*log(xai2)/2.
            enddo
            xa2= exp(sum2/sum1)**2
         endif
      endif
!///////////////
!      write(*,'(a, 1p4g15.7)')
!     *    'xa2 ',  massratio2, beta2, xa2, xa2temp
!c//////////////
      end
!
      subroutine epkaic2(mediax, xc2)
      use modMultipleScat
      implicit none
#include "Zglobalc.h"
!#include "ZepTrackv.h"
#include "Zmedia.h"
      

!
!   note: we neglect atomic electron contribution because it is
!         considered in Moller or Bhabha scattering.
!
!         compute Xc^2 = 4Pi r_0^2 N0 z^2 rho Z^2/A  * integral
!        0 to leng of 1/beta**4/gamma**2  (radian^2)/massratio2
      
      type(epmedia):: mediax    ! input. media
      real*8 xc2   ! output. Xc^2 in radian^2


      real*8 r0, avoganum, const, large
      parameter (r0=2.817d-15, avoganum=6.022d23)
      parameter (const = 4.*pi* r0**2 * avoganum*1.d3, 
     *    large = (pi/2)**2)
!
!
!      integeral 0 to leng of  1/(beta**4 E**2)
!      is approximated as  leng/(beta1**2 gamma1 beta2**2 gamma2)
!      Note: bata**2 *gamma = gamma - 1/gamma
!
      if(gbeta2 .le. 1.d-7) then
         xc2 = large
      else
!     we muse use   rho in kg/m^3  length  in m because of Cosmos
!     but they are cgs here so
!      kg/m^3  m = kg/m^2 =  1000 /10000 g/cm2
!         rho*1000    dl/100.  = rho*dl*10
!c         xc2 = const* cTrack.p.charge**2 * mediax.Z2eff/mediax.Aeff *
!     *         mediax.rho*mediax.rhoc* Move.dl *10./gbeta2 *
!     *         massratio2
!           mediax.Z2eff/mediax.Aeff  = sum(niZi^2)/sum(niAi)
!                we employ Z(Z+1) regardless of e or heavy
!         me^2 ( z/p*beta)^2 = (me/M)^2 z^2/ (g*beta^2)^2
!              < (g*beta^2)^2> --> gbeta2
!     xc2 =const*cTrack%p%charge**2 *  massratio2 *
!     *     (mediax%Z2+mediax%Z)/mediax%A *
!     *     mediax%rho*mediax%rhoc* Move%dl *10./gbeta2
         
         xc2 =const*charge**2 *  massratio2 *
     *        (mediax%Z2+mediax%Z)/mediax%A *
     *        mediax%rho * mediax%rhoc* dlcm *10.0/gbeta2
!     *        mediax%rho * mediax%rhoc* dlcm * 1.0d-5/gbeta2  ; this xxx
!!         xc2 =const*charge**2 *  massratio2 *
!!     *        (mediax%Z2+mediax%Z)/mediax%A *
!!     *        dxgpcm2/gbeta2  --> dxgpcm2*10/gbeta2

         
      endif
      end

      subroutine epElsepaCondense(sgr, theta)
      use modMCScontrol
      use modsoftMCS
      use modcMCS
      implicit none
!     #include "ZepTrackv.h"
#include "ZmediaLoft.h"

      real(8),intent(in):: sgr  ! path length in g/cm2
      real(8),intent(out):: theta ! sampled angle
      
      real(8)::cm2topgrm,s0,s1, s2, ls1, ls2 
      real(8):: mu, cosa, sina
!!!      call ciniFsoft(KEeV, sgr)   is for integration
!!!    from 0 to muc.  We need here from 0 to 1. so
!!!   do below
!            2nd arg. 1 is to get s in cm2  KEeV has been given
!     in modMCScontrol via cfixMCSmodel
!        Cosmos: it is called from cpreSampEIntL.f
!        Epics:  it is called from epMCShard in epMCShard.f
!         which is called from  eppoe in epgen.f. It is placed
!         after cepSampEIntL.
!     TTXSnow  is in modcMCS in cmkMCSconst.f and fixed in
!        cfixMixedConst in cmkMCSconst.f which is in
!        LibLoft/EM/FromCos/MixedMCS/Sampling.
!     cfixMixedConst is called
!     Cosmos:   from cpreSampEIntL. *** in its calling, the
!         first argument must be MediaNo now. (1 could be wrong in V9.00)
!          ==> updated.
!     Epics:   epMCShard       
      
      call cTPXS(TPXSnow, 1,  KEeV, s0, s1, s2)
      cm2topgrm = Media(MediaNo)%mbtoPgrm / 1d-27 
      ls1 =1.0/(s1*cm2topgrm)
      ls2 =1.0/(s2*cm2topgrm)

      avemu = (1.0d0 - exp(-sgr/ls1))/2
      avemu2 = avemu - (1- exp(-sgr/ls2) )/6.d0
      if( avemu > 0.49999999d0) then
!            mu is uniform see bottom
         asoft = 0.5d0
         bsoft = 0.5d0
      else
         bsoft =( 2*avemu - 3*avemu2 )/(1-2*avemu)
         asoft = (1-2*avemu) + bsoft
      endif

      call csampSoftMCSang(mu, cosa)
      theta = acos(cosa) 
      end
!   KEeV-------------=   103204.716982271
!        avemu, avemu2
!      0.499999999066371       0.333333332399704
!              asoft, bsoft
!      0.499999942409857       0.499999940542599
!  ls1, ls2, sgr
! 2.557051000529629E-3  1.121752801379088E-3 5.139364343647575E-2  
!
! KEeV-------------=   96143.8824599196
! avemu, avemu2
! 0.499999999016370       0.333333332349703
! asoft, bsoft
!    0.500000030184747       0.500000028217486
!  ls1, ls2, sgr
! 2.293085720720103E-3 1.012315784158134E-3  4.596862666886253E-02  

      subroutine epPreScatCalc
!     compute g* beta^2,  beta^2, (Me/M)^2
      use modMultipleScat
      implicit none
#include "Zglobalc.h"
#include "Zmass.h"
#include "Zcode.h"
! #include "ZepTrackv.h"

      


      real(8):: g1, g2, g

      if( ready ) return     

      
      g1 = Ebef/mass
      g2 = Eaft/mass
      beta2 = 1.d0 - 1.d0/g1/g2   !        < beta^2>
      gbeta2=(g1 - 1.d0/g1) * (g2 - 1.d0/g2)   ! < (g*beta^2)^2 >
!      gbeta2 =( sqrt(g1*g2)*beta2)**2
!      gbeta2 =g1*g2*(beta2)**2
!      if(gbeta2 == 0.) then
!         g = (g1+ g2)/2.
!         gbeta2 = (g-1.d0/g)**2
!      endif

      if( code  /= kelec ) then
         massratio2 = (masele/mass)**2 ! (me/m)^2
      else
         massratio2 = 1.0d0
      endif

      ready = .true.
      end


