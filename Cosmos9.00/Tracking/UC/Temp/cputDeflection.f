      subroutine cputDeflection
       implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"


       real*8 dt, dispx, dispy
       logical nodef
       if(Reverse .eq. 0) then
          nodef = OneDim .ne. 0  .or.
     *    (mod(HowGeomag, 10) .eq. 1 .and. Zfirst%pos%depth .eq. 0.)
       else
          nodef = .false.
       endif
       if(.not. nodef) then
!
!         compute   magnetic deflection first. independently of scattering.
!         system is xyz;  if Efield exists, together with it
          call cmagDef
       endif
       if(.not. nodef .and. Reverse .eq. 0) then
!                this is for multiple scattering
          call celecDef(dispx, dispy)  ! effect alrady put in
                    ! MovedTrack. dispx, y are dummy
       endif
       if(.not. nodef) then
          call csetPos(MovedTrack%pos)
          call cgetZenith(MovedTrack, MovedTrack%vec%coszenith)
!          reset momentum to be compatible with direction cos.
          call cresetMom(MovedTrack)
!
          if(TimeStructure .and. Reverse .eq. 0) then  ! only for mul. scat
!                compute excess path length to be added to streight path
             call cexcessLen(dispx, dispy, dt)
             if(Beta .gt. 0.) then
                MovedTrack%t = MovedTrack%t + dt/Beta
             endif
          endif
       endif
       end
!      *************************
#include "ZsubstRec.h"
       subroutine cmagDef
       use modEfield
       use modSetIntInf
       implicit none
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"

!

!              by  Geomag  (dr and ddirec)
       type(coord)::dispmr, dispmd
       type(magfield)::tempmag
       real*8  leng
       real*8 temp1, temp2
!cc    *   ,norm
       integer i

       integer icon
       type(coord):: middle
       logical usemiddle, high, usecmiddle
       real*8 cheight
       data cheight/30.d3/  
       real*8 newmom(3), dE, E1, E2, p1sq, p2sq, newE, norm
       save cheight

       leng = IntInfArray(ProcessNo)%length
       if( HowEfield >= 1 ) then
          call cdefByMagAndE(TrackBefMove,  leng,  dispmr, dispmd,
     *     newmom)
          p1sq =dot_product(TrackBefMove%p%fm%p(1:3),
     *           TrackBefMove%p%fm%p(1:3))
          E1 = sqrt( p1sq + TrackBefMove%p%mass**2)
          p2sq =dot_product(newmom(1:3),newmom(1:3))
          E2=sqrt(p2sq +  TrackBefMove%p%mass**2)
          dE=E2-E1  !  may be + or -  depending on Ef and charge
          ! dE/dx loss has been put in MovedTrack. adjust it again
          
          newE = MovedTrack%p%fm%p(4) + dE
          if( newE <  TrackBefMove%p%mass ) then
             newE =  TrackBefMove%p%mass
             newmom(:) = 0
          else
             ! keeping the direction,adjust new momentum to be
             ! consistent with newE
             if( p2sq <= 0. ) then
                write(0,*) ' p2sq = ',p2sq
                write(0,*)
     *          ' TrackBefMove%p=',TrackBefMove%p%fm%p(1:4)
                write(0,*) ' dE, newE=',dE, newE
                write(0,*) ' leng =', leng
                write(0,*) 'dispmr=',dispmr%r(:)
                write(0,*) 'dispmd=',dispmd%r(:)
                write(0,*)
     *            'code, chg=',TrackBefMove%p%code,
     *            TrackBefMove%p%charge

                newE = TrackBefMove%p%mass
                newmom(:) = 0
                dispmd%r(:) = (/0.,0.,1./)
             else
                norm = sqrt( (newE**2 - TrackBefMove%p%mass**2)/p2sq)
                newmom(:) = newmom(:) *norm
             endif
          endif
          MovedTrack%p%fm%p(4) = newE
!          note. dispmr is r(new)-r(old) vector and different
!          from other routines below. 
!           dispmd is new dir and set at 100
          MovedTrack%pos%xyz%r(:) =  TrackBefMove%pos%xyz%r(:) +
     *         dispmr%r(:) 
          MovedTrack%p%fm%p(1:3) = newmom(:)
          goto 100
       endif
!
!            UseRungeKutta         Height       
!                0                  any         Mag and cmagneticDef

!                1                 >cheight     middle Mag and cmagneticDef
!                1                 <            Mag and cmagneticDef

!                2                 >            middle Mag and cmagneticDef 
!                                               or cbDefByRK2        
!                2                 <            Mag and cmagneticDef 

!                3                 >            middle Mag and cmagneticDef 
!                                               or cbDefByRK        
!                3                 <            Mag and cmagneticDef 

!                4                 >            use Mag at curved middle point
!                                               and estimate end point by 
!                                               cmagneticDef or cbDefByRK2
!                4                 <            Mag and cmagneticDef   

!                5                 >            use Mag at curved middle point
!                                               and estimate end point by 
!                                               cmagneticDef or cbDefByRK
!                5                 <            Mag and cmagneticDef   

!                6                 >            cbDefByRK2
!                6                 <            Mag and cmagneticDef 

!                7                 >            cbDefByRK
!                7                 <            Mag and cmagneticDef 
!
!                8      at any height           cbDefUser ; interface is
!                                               same as cbDefByRK
!
       if(UseRungeKutta .eq. 8 ) then
          call cbDefUser(TrackBefMove,  leng,  dispmr, dispmd,
     *    MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
          goto 100
       endif
       high = TrackBefMove%pos%height .gt. cheight
       usemiddle = ( UseRungeKutta .ge. 1 .and.
     *               UseRungeKutta .le. 3 .and.
     *               high )
       usecmiddle = ( UseRungeKutta .ge. 4 .and.
     *               UseRungeKutta .le. 5 .and.
     *               high )


       if(usemiddle) then
          do i = 1, 3
             middle%r(i) = TrackBefMove%pos%xyz%r(i)+ 
     *          leng/2 * TrackBefMove%vec%w%r(i)
          enddo
          middle%sys='xyz'

          call cgeomag(YearOfGeomag,  middle,
     *                tempmag, icon)
          call ctransMagTo('xyz', middle,
     *        tempmag, tempmag)
       elseif(usecmiddle) then
!            get curved middle point
          call cmagneticDef(TrackBefMove, Mag, leng/2.0d0,
     *     dispmr, dispmd)  
          do i = 1, 3
             middle%r(i) = TrackBefMove%pos%xyz%r(i) + dispmr%r(i)
          enddo
          middle%sys ='xyz'
          call cgeomag(YearOfGeomag,  middle,
     *                tempmag, icon)
          call ctransMagTo('xyz', middle,
     *        tempmag, tempmag)
       endif

       if(usemiddle .and. UseRungeKutta .eq. 1) then
          Mag = tempmag
       elseif( usemiddle .and.  UseRungeKutta .eq. 2 ) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.4 %, use RungeKutta.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.0d-3) then
             call cbDefByRK2(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       elseif(usemiddle .and. UseRungeKutta .eq. 3) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.45 %, use RungeKutta.
!                  next 4.5d-3 is very sensitive to cpu time
!                 if it was 4.0d-3, cpu time becomes twice.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.5d-3) then
             call cbDefByRK(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       elseif(usecmiddle .and. UseRungeKutta .eq. 4) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.1 %, use RungeKutta.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.0d-3) then
             call cbDefByRK2(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       elseif(usecmiddle .and. UseRungeKutta .eq. 5) then
          temp1 = tempmag%x**2 + tempmag%y**2 + tempmag%z**2
          temp2 = Mag%x**2 + Mag%y**2 + Mag%z**2
!           if  dB/B~ dB^2/B^2/2 > 0.45 %, use RungeKutta.
          if( abs(temp1-temp2)/temp1/2.0 .gt. 4.5d-3) then
             call cbDefByRK(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC             
             MovedTrack%pos%xyz = dispmr ! this is not a dispalcement vector
#else 
             call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
             goto 100
          else
             Mag = tempmag
          endif
       endif
!!!!!!   default comes here
       if(UseRungeKutta .le. 5 .or. .not. high ) then
          call cmagneticDef(TrackBefMove, Mag, leng,
     *     dispmr, dispmd)

          do i = 1,  3
             MovedTrack%pos%xyz%r(i) = TrackBefMove%pos%xyz%r(i) +
     *         dispmr%r(i)
          enddo
       elseif(UseRungeKutta .eq. 6) then
          call cbDefByRK2(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC
          MovedTrack%pos%xyz = dispmr  ! this is not a dispalcement vector
#else 
          call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
                                      !  but the vector itself.
       else
          call cbDefByRK(TrackBefMove,  leng, dispmr, dispmd,
     *         MovedTrack%p%fm%p(1:3))
#if defined SUBSTREC
          MovedTrack%pos%xyz = dispmr  ! this is not a dispalcement vector
                                      !  but the vector itself.
#else 
          call csubstcoord(dispmr, MovedTrack%pos%xyz)
#endif
       endif
 100   continue
#if defined SUBSTREC
       MovedTrack%vec%w = dispmd
#else
       call csubstcoord(dispmd, MovedTrack%vec%w)
#endif
       end
!      ==============================================================
      subroutine celecDef(dx, dy)
      use modSetIntInf
      implicit none

#include  "Ztrack.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"
!!!!!!!!!!!
       real(8),save:: X0=365.667
!!!!!!       
      real*8 dx, dy             ! output. scatttering displacement
!     
      type(coord)::dircos
!              by  Multiple Scattering
      type(coord)::dsa   ! dire ccos of scattering angle
      type(coord)::w
      real*8  sint, cs, sn, tmp, avx, avy, disp, dt, dl
      real*8 r, g1, g2, gf1, gf2, beta2, tetarms
      real*8 theta
      
      call cmulScat(theta)
      if(theta .lt. 0.01d0) then
!                 cos
!     dsa%z = 1.-theta**2/2
         dsa%r(3) = 1.-theta**2/2
         sint = theta
      else
!     dsa%z = cos(theta)
         dsa%r(3) = cos(theta)
         sint = sin(theta)
      endif
!        azimuthal angle
      call kcossn(cs, sn)
!     dsa%x = sint * cs
!     dsa%y = sint * sn
!     dsa%z = cos(theta)
      dsa%r(1) = sint * cs
      dsa%r(2) = sint * sn
      dsa%r(3) = cos(theta)    
      dt = IntInfArray(ProcessNo)%thickness/ X0 ! r%l
      dl = IntInfArray(ProcessNo)%length !  m

!      if(Moliere .ge. 0) then
!     if( ALateCor ) then       ! >v7.645
!         v>=7.655      
      if( (Moliere == 0 .and. AlateCor>=1) .or. 
     *     AlateCor == 2) then
!     ALateCor D=1.  if Gaussian (Moliere=0), correlation is taken
!     into account.
!     If 0, no correlation is considered.
!        2, correlation based on Gaussian formula is forced even if
!           Moliere is not 0.
         
!         sample displacement correlated to theta
!                 this is the same as P.D.B though look like
!                diff.
        tmp = dl/2.d0
        avx = tmp * dsa%r(1)
        avy = tmp * dsa%r(2)
!                dispersion
        gf1 = TrackBefMove%p%fm%p(4)/MovedTrack%p%mass
        gf2 = MovedTrack%p%fm%p(4)/MovedTrack%p%mass
        beta2 = 1.d0 - 1.d0/gf1/gf2
        if(beta2 .le. 0.) then
           disp = 0.d0
        else
           if(dt .gt. 1.d-3) then
              tetarms = Es/TrackBefMove%p%fm%p(4)*
     *            abs(MovedTrack%p%charge)*
     *            sqrt(dt)*(1.0 + 0.038*log(dt))
           else
              tetarms = Es/TrackBefMove%p%fm%p(4)*
     *                  abs(MovedTrack%p%charge)*
     *              sqrt(dt)
           endif
           disp=tetarms/sqrt(6.d0*beta2)*dl/2.d0
!               sample 2 independent gaussian variables
!             with mean 0 and var 1
        endif
        call kgauss2(0.d0, 1.0d0, g1, g2)
        dx = g1 * disp + avx
        dy = g2 * disp + avy
!                  displacement
        r=sqrt(dx*dx+dy*dy)     ! in m
!              direction cos of vector r in original sys.
        if(r .ne. 0.) then
!           w%x = dx/r
!           w%y = dy/r
!           w%z = 0.
           w%r(1) = dx/r
           w%r(2) = dy/r
           w%r(3) = 0.d0
           
!                 transform wx,wy,wz to original sys.
!                    TrackBefMove is better
           call ctransVectZ(TrackBefMove%vec%w, w, w)
!               r is already in m.
!              add scattering effect.
!              r*w is displacement by scattering
!           MovedTrack%pos%xyz%x = r*w%x + MovedTrack%pos%xyz%x
!           MovedTrack%pos%xyz%y = r*w%y + MovedTrack%pos%xyz%y
!           MovedTrack%pos%xyz%z = r*w%z + MovedTrack%pos%xyz%z
           MovedTrack%pos%xyz%r(1:3) =
     *         r*w%r(1:3)+ MovedTrack%pos%xyz%r(1:3)
        endif
      else
         dx = 0.
         dy = 0.
      endif
!        convert scattering angle at end of path to
!        original system . MovedTrack is better since
!        mag. def is contained there already.
      call ctransVectZ(MovedTrack%vec%w, dsa,
     *      MovedTrack%vec%w)

      end
!     ****************************************
!     *                                                          *
!     * cmulScat: multiple Coulomb scattering
!     *                                                          *
!     ****************************************
!
!
      subroutine cmulScat(theta)
!      use modXsecMedia
      use modMCScontrol
!!!!!!      
      use modSetIntInf
!!!!!!!!!!      
      implicit none
!       Using  cTrack and Move, compute scattering angle 
!
!!!!!!!!!11
#include "ZmediaLoft.h"
!!!!!!!!!      
#include "Zglobalc.h"
#include  "Zcode.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"

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
            call cSampMol(Media(1), theta, cond)  ! rigorous Moliere
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
      subroutine cmulScat2(theta)
      use modSetIntInf
      implicit none
!       Gaussian scattring
!
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
#include  "Zelemagp.h"

      

      real*8 theta ! output. sampled spatial angle in radain.

      real*8 tetarms, g1, g2, u, beta2, dt
      integer nc
      
      real*8 hpi 
      parameter(hpi = pi/2.)

!!!!!!1      
      real(8),save:: X0=365.667
!!!!!!1      
      g1 = TrackBefMove%p%fm%p(4)/MovedTrack%p%mass
      g2 = MovedTrack%p%fm%p(4)/MovedTrack%p%mass
      beta2 = 1.d0 - 1.d0/g1/g2
      if(beta2 .le. 0.) then
         tetarms = 0.
      else
         dt = IntInfArray(ProcessNo)%thickness/ X0  ! r%l
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
      use modMCScontrol
      use modsoftMCS
      use modcMCS
      use modSetIntInf
!      use modXsecMedia
      implicit none
!!!!!!!!!
#include "ZmediaLoft.h"
!!!!!!!      
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
      cm2topgrm = Media(1)%mbtoPkgrm/ 1d-27 
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
      subroutine cmoliere(rhoin, 
     *           z, mass, g1, g2,  leng, teta, cond)
      implicit none
!        Moliere theory of multiple scattering angle.
!
      real*8 rhoin  ! input. average media density in kg/m^3

      integer z  ! input.  charge of the particle is ze
      real*8 mass ! input.  mass of the particle in GeV
      real*8 g1  ! input.  gamma factor at the path head
      real*8 g2  ! input.  gamma factor at the path end
      real*8 leng ! input. length the charged particle travelled in m.      
      real*8 teta ! output. sampled spatial angle in radain.
      integer cond ! output. 0 ok. non-0. Moliere theory not applicable

      real*8  rho, gbeta2, beta2, massratio2
      common /Zcmedia/ rho, gbeta2, beta2,
     *       massratio2

!     *********************
      real*8 xc2, xa2, bp, b, u
      real*8 a0, a1, a2,  sum, ra2, ra2inv

      integer icon
      real*8 rejf1, rejf21, rejf22
      real*8 x
!       
!      rejection function for redueced angle < 1.8
!          x is suqre of reduced angle better than 0.2 %
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
      rho = rhoin
!        bgeta2   (g-1/g) = g*beet^2 = g*(1-1/g^2)
      gbeta2=(g1 - 1.d0/g1) * (g2 - 1.d0/g2)
      beta2 = 1.d0 - 1.d0/g1/g2
      massratio2= (0.511d-3/mass)**2  ! (me/m)^2
!    ............................

!          get Xc^2
      call ckaic2(z,  leng,  xc2)
!          get Xa^2
      call ckaia2(z,  xa2)

!          b -log(b) = b'
      bp = log(xc2/xa2/1.167)

      if(bp .lt. 3.395) then
!         Moliere theory cannot be appliled; use Gaussian later (almost no scattering)
         cond = 1
      else
         cond = 0
         call cblogb(bp, b, icon)
         a0 = max(1.d0 - 5/b, 0.d0)  ! use single scattering term if b<=5. 
         icon = 1               ! make 0 if no rejection
!                the sampling function decomposition is explained in Test/....tex
         do while (icon .ne. 0)
            a1 = 5.21062/b
            a2 = 0.7128/b
            sum = a0 + a1 + a2
            call rndc(u)
            if(a0/sum  .gt. u) then
!             sample reduced angle from exp(-x) dx where x = reduced
!               angle^2.
               call rndc(u)
               ra2 = -log(u)
               icon =0
            elseif( (a0+a1)/sum .gt. u) then
!            sample reduced angle from exp(-x) dx (same as above but
!                  in the region of ra < 1.8
               call rndc(u)
               ra2 = -log(1.-u/1.04076)
!                rejection function
               call rndc(u)
               if(u .lt. rejf1(ra2)) then
                  icon = 0
               endif
            else
!             sample reduced angle from 2xc2 x^-4dx
               call rndc(u)
               ra2 = 3.24/u
!               rejection function
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
         enddo
         teta =sqrt( ra2 * xc2 * b)
      endif
      end
!     ***************************
      subroutine cblogb(c, b, cond)
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
      endif
      end
      subroutine ckaia2(z,   xa2)
      implicit none
#include "Zair.h"
!        compute Xa^2; assume the Xa^2 is weakly
!       dependent on Z, we use average Z=zave for
!       calculation.
!
      integer z  ! input. charge of the charged particle is ze
      real*8 xa2   !  output. Xa^2.

      real*8  rho,  gbeta2, beta2, massratio2
      common /Zcmedia/ rho, gbeta2, beta2, 
     *       massratio2
      
      real*8   alpha, const, pi, large
      parameter (alpha = 1./137., const = (1.13*alpha)**2 )
      parameter (pi = 3.1415, large = (pi/2.)**2)

      if(gbeta2 .le. 0.) then
         xa2 = large
      else
!              since air Z is close N and O's Z, we
!           use simply average of Z here.
         xa2 = const * TargetZ2_3rd * massratio2 *
!     *   (1.13 * beta2 + 3.76*(alpha*z*zave)**2) /gbeta2
     *   (1.13 *beta2 + 3.76*(alpha*z*TargetAtomicN)**2)
     *     /gbeta2
      endif
      end
!
      subroutine ckaic2(z, leng, xc2)
      implicit none
#include "Zair.h"
!
!   note: we neglect atomic electron contribution because it is
!         considered in Moller or Bhabha scattering.
!
!         compute Xc^2 = 4Pi r_0^2 N0 z^2 rho Z^2/A  * integral
!        0 to leng of 1/beta**4/gamma**2  (radian^2)/massratio2
!
      integer z     ! input. charged particle charge is ze.
      real*8 leng  ! input. length traveled by the charge particle in m
      real*8 xc2    ! output. Xc^2 in radian^2

      real*8  rho, gbeta2, beta2, massratio2
      common /Zcmedia/ rho, gbeta2, beta2,
     *       massratio2

      real*8 r0, avoganum, const, large, pi
      parameter (r0=2.817d-15, avoganum=6.022d23, pi= 3.1415)
      parameter (const = 4.*pi* r0**2 * avoganum*1.d3, 
     *    large = (pi/2)**2)
!
!
!      integeral 0 to leng of  1/(beta**4 E**2)
!      is approximated as  leng/(beta1**2 gamma1 beta2**2 gamma2)
!      Note: bata**2 *gamma = gamma - 1/gamma
!

      if(gbeta2 .le. 0.) then
         xc2 = large
      else

         xc2 = const* z* z * massratio2 *
!                    < Z(Z+1)> /<A>
     *   (TargetZ2 + TargetAtomicN)/TargetMassN 
     *   * rho * leng/gbeta2   
      endif
      end
!        This a copy of epSampMol.f. MOdificaiton is
!           epxxx is changed to cxxx
!           cTrack --> TrackBefMove
!           Move --> MovedTrack
!           Move.dx-->IntInfArray(ProcessNo).thickness/10.
      subroutine cSampMol(mediax, teta, cond)
!      use modXsecMedia
      
      implicit none
#include "Zmedia.h"      
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
!#include "Zglobalc.h"
!#include "Zmass.h"
!#include "Ztrack.h"
!#include "Ztrackv.h"
      
!      type(epmedia):: mediax  ! input. media 
      type(epmedia),intent(in):: mediax 
      real(8),intent(out)::teta !  sampled spatial angle in radian.
      integer,intent(out)::cond !  0 ok. non-0. Moliere theory not applicable
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
         call cSampMol0(mediax, teta, icon)
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
      
      
      subroutine cSampMol0(mediax, teta, cond)
!      use modXsecMedia
      use cmodSampMolReducedA
      implicit none
#include "Zmedia.h"
!        Moliere theory of multiple scattering angle.
!        with improvement by Bethe.
!#include "Zglobalc.h"
!#include "Zmass.h"
!#include "Ztrack.h"
!#include "Ztrackv.h"

!      type(epmedia):: mediax  ! input. media 
      type(epmedia),intent(in):: mediax

      real(8),intent(out)::teta !  sampled spatial angle in radian.
      integer,intent(out)::cond !  0 ok.  Moliere used
                                !  1 ok.  G.S is used.
                                !  2 n.g  Gauss should be used.

      real*8  B, xc2, Bxc2
      common /Zmoliere/ B, xc2, Bxc2

      real*8  gbeta2, beta2, massratio2
      common /Zcmedia/  gbeta2, beta2,  massratio2

      
!     *********************
      real*8  expb
      real(8):: sb,  tetasq

      call cPreScatCalc  ! get g*beta^2, beta^2, (Me/M)^2 in Zcmdeia
!          get Xc^2
      call cXaic2(mediax,   xc2)
!           we don't get xa2 but directly get exp(b)
      call cMolExpb(mediax,   expb)
      sb =log(expb)  ! 
!///////////
!      write(*,'(a, 1p,6g14.5)')
!     *  ' sb ', sb, beta2, Move.dt, Move.dl, cTrack.p.fm.p(4),
!     *          Move.track.p.fm.p(4)
!//////////////
!          B -log(B) = b  
!///      if(sb .lt. 3.395d0) then
      if(sb .lt. 1.5d0) then
!         Moliere theory cannot be appliled; use Gaussian later
!         (almost no scattering)
         cond = 2
      else
         call cMolBlogB(sb, B)  !  B
         Bxc2 = B*xc2
         if(Bxc2 > 1.0d0) then
!             normally this happens for low E, which stops soon
!             B  so not serious.
!               use GS method ? --> cond=1
            cond =2   ! use Gauss at present
         else
            call cSampMolReducedA(B, tetasq, cond)
            if( cond == 0 ) then   ! now alwasy cond=0
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
      subroutine cXaic2(mediax, xc2)
!       compute xc2 of Moliere's theory of scattring.
!      epPreScatCalc must have been called beforehand for 
!      each scattering.
!     use modXsecMedia
      use modSetIntInf
      implicit none
#include "Zglobalc.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zcode.h"
#include "Zmedia.h"
      
!      type(epmedia):: mediax ! input. media
      type(epmedia),intent(in):: mediax

      real*8 xc2   ! output. Xc^2 in radian^2

      real*8  gbeta2, beta2, massratio2
      common /Zcmedia/ gbeta2, beta2, massratio2


!
!      integeral 0 to leng of  1/(beta**4 E**2)
!      is approximated as  leng/(beta1**2 gamma1 beta2**2 gamma2)
!      Note: bata**2 *gamma = gamma - 1/gamma
!


!      if( cTrack.p.code == kelec) then
      if( TrackBefMove%p%code == kelec) then
         xc2 = mediax%MoliereForXc2
      else
!         xc2 = 0.6011d0*cTrack.p.charge**2 * massratio2 *
         xc2 = 0.6011d0*TrackBefMove%p%charge**2 * massratio2 *
     *     mediax%Z2/mediax%A 
      endif
!         Move.dx is for EPICS. 
!      xc2 = xc2*Move.dx/gbeta2  ! dx is in g/cm2. no need to consider rhoc.
                                ! real length is dx/(rho*rhoc)
!       thickness is kg/m2;*  1000/10000 = 1/10  g/cm2
      xc2 = xc2*IntInfArray(ProcessNo)%thickness/10./gbeta2  ! 

      end

      subroutine cMolExpb(mediax, expb)
!            compute exp(b) of Moliere's scattering theory.
!            epPreScatCalc must have been called beforehand for 
!     each scattering.
      use modSetIntInf 
!      use modXsecMedia
      implicit none
#include "Zglobalc.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zmedia.h"
      
!      type(epmedia):: mediax ! input. media
      type(epmedia),intent(in):: mediax
      real*8 expb   ! output. exp(b) of Moliere's theory.

      real*8  gbeta2, beta2, massratio2
      common /Zcmedia/ gbeta2, beta2, massratio2
      
      real(8):: sum, alfai2, temp
      integer::i 

      if( beta2 > 0.99 ) then  ! beta > 0.995 --> regards 1
         sum = mediax%MoliereExpb
      else
         sum =0.
!         temp = (cTrack.p.charge/137.0)**2/beta2
         temp = (TrackBefMove%p%charge/137.0)**2/beta2
         do i = 1,  mediax%noOfElem

            alfai2 = mediax%elem(i)%Z**2 *temp

            sum = sum +
     *       mediax%No(i) *
     *       mediax%elem(i)%Z**0.3333*(mediax%elem(i)%Z + 1.d0)
     *       / (1.d0 + 3.327d0*alfai2)
         enddo
         sum = 6702.d0 *sum/mediax%A
      endif
!      expb = sum *Move.dx * cTrack.p.charge**2 /beta2
      expb = sum *IntInfArray(ProcessNo)%thickness/10. *
     *   TrackBefMove%p%charge**2 /beta2
     
      end


      subroutine cPreScatCalc
!           compute g* beta^2,  beta^2, (Me/M)^2
      implicit none
#include "Zglobalc.h"
#include "Zmass.h"
#include "Zcode.h"
#include "Ztrack.h"
#include "Ztrackv.h"

      real*8  gbeta2, beta2, massratio2
      common /Zcmedia/  gbeta2, beta2,  massratio2

      real(8):: g1, g2, g

!      g1 = cTrack.p.fm.p(4)/cTrack.p.mass
      g1 = TrackBefMove%p%fm%p(4)/TrackBefMove%p%mass
!      g2 = Move.Track.p.fm.p(4)/cTrack.p.mass
      g2 = MovedTrack%p%fm%p(4)/TrackBefMove%p%mass
      beta2 = 1.d0 - 1.d0/g1/g2   !        < beta^2>
!      gbeta2=(g1 - 1.d0/g1) * (g2 - 1.d0/g2)   ! < (g*beta^2)^2 >
!      gbeta2 =( sqrt(g1*g2)*beta2)**2
      gbeta2 =g1*g2*(beta2)**2
!      if(gbeta2 == 0.) then
!         g = (g1+ g2)/2.
!         gbeta2 = (g-1.d0/g)**2
!      endif

!      if(cTrack.p.code == kelec) then
      if(TrackBefMove%p%code == kelec) then
         massratio2= 1.
      else
!         massratio2= (masele/cTrack.p.mass)**2 ! (me/m)^2
         massratio2= (masele/TrackBefMove%p%mass)**2 ! (me/m)^2
      endif
      end
