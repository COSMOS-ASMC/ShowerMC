      subroutine epdeflection(icon)
      use moddebug
      use modMCScontrol
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"
!        All deflections are managed here:
!        Move.Track is assumed to be the streight extension 
!        of cTrack.  Here deviation from the initial Moved.Track
!        must be computed and  added to it successively.
!        However, if we compute magnetic deflection, we can
!        compute the end point coordinate and direction cosines
!        directly.  So we compute Magnetic def. first and
!        add other effect later.
!
!     1) magnetic deflection
!     2) electric deflection
!     3) multiple Coulomb scatteringc
      integer icon  ! output.  icon=1 when infinite loop or unchanged

       type(epmove):: moveSave
      integer nx

!          save current Move
      moveSave%Track%pos = Move%Track%pos
      if(MagField .ne. 0) then
         if(Bfield%x .ne. 0. or. Bfield%y .ne. 0. .or.
     *      Bfield%z .ne. 0. ) then
            call epdefByB(icon)
         endif
      endif
      if(ElecField .ne. 0) then
         call epdefByE
      endif
      if(Media(MediaNo)%name .ne. 'sp' .and.
     *     Media(MediaNo)%name .ne. 'hollow') then
         if(( Move%Track%p%fm%p(4)- Move%Track%p%mass > 
     *       MaxErgMCS) .or. ActiveMCS == "El_con" .or.
     *        ActiveMCS == "Mol") then
            call epdefByMScat
         endif
      endif
!       check if position moved to another component
!       due to scattering, if so, resotore the position
!       without  scattering/deflection
      if(.not. Move%Cross) then  ! ++++++++++
         call eppos2cn(Cn, Move%Track, nx)
         if(Cn .ne. nx) then
!           the new pos. dose not remain in the same component
!            no displacement due to scattering
            Move%Track%pos = moveSave%Track%pos
         endif
      endif
      end
!     ***********************
      subroutine epdefByMScat
!         deflection by Coulomb scattering
      implicit none
#include "ZepTrackv.h"
#include "ZepTrackp.h"

       type(epDirec)::  dsa   ! dire ccos of scattering angle
       type(epDirec)::  w
      real*8 theta, sint, cs, sn, tmp, avx, avy, disp
      real*8 r, dx, dy, g1, g2, gf1, gf2, beta2, tetarms
      real(8):: sinx, cosx


      call epmulScat(cTrack%p, Move%Track%p, Move%dl, Move%dx,
     *      Media(MediaNo), theta)    ! theta may > pi/2 for Moliere=T or >=2

      if(theta .lt. 0.3d0) then
!                 cos
         dsa%z = 1.-theta**2/2   ! 0.04 % error
         sint = (-theta**2/6.d0 +1.d0)*theta 
      else
         dsa%z = cos(theta)
         sint = sin(theta)
      endif
!        azimuthal angle
      call kcossn(cs, sn)
!
      dsa%x = sint * cs
      dsa%y = sint * sn

      if( Move%Cross ) then
!        dont consider ang. late correl.
      elseif( ( Moliere==0 .and. ALateCor >= 1) .or.
     *        ALateCor == 2 ) then
!         sample displacement correlated to theta 
!                 this is the same as P.D.B though look like
!               diff.
         tmp = Move%dl/2.d0
         avx = tmp * dsa%x
         avy = tmp * dsa%y
!                dispersion
         gf1 = cTrack%p%fm%p(4)/cTrack%p%mass
         gf2 = Move%Track%p%fm%p(4)/cTrack%p%mass
         beta2 = 1.d0 - 1.d0/gf1/gf2
         if(beta2 .le. 0.) then
            disp = 0.d0
         else
            if(Move%dt .gt. 1.d-3) then
!              beta is considered later in disp
               tetarms = Escat/cTrack%p%fm%p(4)*abs(cTrack%p%charge)*
     *      sqrt(Move%dt)*(1.0 + 0.038*log(Move%dt))
!     *      sqrt(Move.dt)

            else
               tetarms = Escat/cTrack%p%fm%p(4)*abs(cTrack%p%charge)*
     *              sqrt(Move%dt)
            endif
            disp=tetarms/sqrt(6.d0*beta2)*Move%dl/2.d0
!               sample 2 independent gaussian variables
!             with mean 0 and var 1
         endif
         call kgauss2(0.d0, 1.0d0, g1, g2)
         dx = g1 * disp + avx
         dy = g2 * disp + avy

!                  displacement

         r=sqrt(dx*dx+dy*dy)

!              direction cos of vector r
         if(r .ne. 0.) then

            sinx = r/Move%dl
            if( Zcorrec  .and.  sinx < 1.0d0) then
 !             In the system of which z is the same as
 !              current particle Z direction,
 !             current ptcl pos=(0,0,Move.dl)
!        Due to the displacement, the length, Move.dl,  is
!            unchaged but the, x,y,z must be modified.
!          to (dx, dy, Move.dl*cosx)  where cosx is comes from
!          sinx= r/Movd.dl
               if(sinx < 0.1) then               
                  cosx = 1.0d0- sinx**2/2
               else
                  cosx = sqrt(1.d0-sinx**2) 
               endif

               w%x = dx
               w%y = dy
               w%z = Move%dl*cosx
                  !  transform w to the original system
               call eptransVect(cTrack%w, w, w)
                                     ! not Move.Track: c.f Zcorrec=f
                                     ! case below 
               Move%Track%pos%x = cTrack%pos%x +
     *           w%x*Move%dl
               Move%Track%pos%y = cTrack%pos%y +
     *           w%y*Move%dl
               Move%Track%pos%z = cTrack%pos%z +
     *           w%z*Move%dl
            else
               w%x = dx/r
               w%y = dy/r
               w%z = 0.
!                 transform w1,w2,w3 to original sys.
               call eptransVect(cTrack%w, w, w)
!               r is already in cm.
!              add scattering effect.
!              r*w is displacement by scattering
                                      ! this is Move
               Move%Track%pos%x = r*w%x + Move%Track%pos%x
               Move%Track%pos%y = r*w%y + Move%Track%pos%y
               Move%Track%pos%z = r*w%z + Move%Track%pos%z
            endif
         endif
      endif
!        convert scattering angle at end of path to
!        original system
      call eptransVect(Move%Track%w, dsa, Move%Track%w)
      end
!     ********************
      subroutine epdefByB(icon)
#include "ZepTrackp.h"
#include "ZepTrackv.h"
       integer icon  ! output. =1 if infinite loop or unchaged
       type(ep3Vec)::  dispmr
       type(epDirec)::  dispmd
       real*8  norm
       real*8 path, radius, loop
       data   path/0./
       save   path, loop
!
!           trick to avoid semi-infinite loop in vaccum 
!
       if(Media(MediaNo)%name .ne. 'sp' .and.
     *    Media(MediaNo)%name .ne. 'hollow'  ) then
          pathInB = 0.
       endif
       if(pathInB .eq. 0.) then
          call epmagDefR(cTrack, Bfield, radius) ! rough radius
!                 if  a ptcl runs  SyncLoop times of cyclotron  radius
!                 discard it
          loop = SyncLoop * radius
       endif
       pathInB = pathInB + Move%dl
       if(pathInB .gt. loop) then
          Move%Trunc = .True.
          icon =1
          pathInB = 0.
       endif

       call epmagneticDef(cTrack, Bfield, Move%dl, dispmr, dispmd)
!
!         dispmr is
!         displacement to be added to cTrack. (not Move.Track)
!              ++++++++++++++
       if(.not. Move%Cross) then
          Move%Track%pos%x = cTrack%pos%x + dispmr%x
          Move%Track%pos%y = cTrack%pos%y + dispmr%y
          Move%Track%pos%z = cTrack%pos%z + dispmr%z
       endif
       Move%Track%w%x = dispmd%x
       Move%Track%w%y = dispmd%y
       Move%Track%w%z = dispmd%z
      end
      subroutine epdefByE
      end

!      Compute magnetic deflection (angle and displacement)
!      
!     ************************************************  
      subroutine epmagneticDef(aTrack, B, leng, dr,  dir)
!     ************************************************  
      implicit none
#include "ZepTrack.h"
#include "Zep3Vec.h"

       type(epTrack)::  aTrack     ! input. a charged ptcl at initial pos.
       type(ep3Vec)::  B            ! input. magnetic field vector. in T
                                  !  in local coordinate
      real*8  leng         ! input. length  travelled in cm
       type(ep3Vec)::   dr  ! output. displacement vector in cm
       type(epDirec):: dir  ! output.  new direction cos 
!    This is exact if B is const.
!   Note that:
!      dr must be added  as  newpos = oldpos +  dr
!      dir is the new direction cos.  These two are differenct
!      from the old version.
!

!
       type(epPos)::  Tx, Ty, Tz, w, r

      real*8 Babs, sint, cost, rgt, sinl, cosl
      real*8   pabs, temp
!
!       get abs(B)
      Babs = sqrt( B%x**2 + B%y**2 + B%z**2 )
      if(Babs .gt. 0.) then
!            p in GeV
         pabs = sqrt(aTrack%p%fm%p(4)**2 - aTrack%p%mass**2)
         if(pabs .gt. 0.) then
!              gyroradius in cm. B in T, p in GeV
!              For negative charge, rgt < 0. This is o.k
!                 real radius is rgt*sintc

            rgt = 3.3358d2*pabs/aTrack%p%charge/Babs

            Tz%x = B%x/Babs
            Tz%y = B%y/Babs
            Tz%z = B%z/Babs

            call cscalerProd(aTrack%w, B, cost)
            cost =max( min(cost/Babs, 1.0d0), -1.d0)
            sint = sqrt(1.0d0-cost**2)
            if(abs(sint) .gt. 1.d-6) then            
               call cvecProd(aTrack%w, Tz, Tx)
               Tx%x = Tx%x/sint
               Tx%y = Tx%y/sint
               Tx%z = Tx%z/sint
               call cvecProd(Tz, Tx, Ty)
!////////////////
!               call checkuv('Tx', sint, Tx)
!               call checkuv('Ty', sint, Ty)
!               call checkuv('Tz', sint, Tz)
!               call checkorth('TxTy', sint, Tx, Ty)
!               call checkorth('TxTz', sint, Tx, Tz)
!//////////////////////               
            else
               Tx%x = 1.d0
               Tx%y = 0.
               Tx%z = 0.
               Ty%x = 0.
               Ty%y = 1.d0
               Ty%z = 0.
            endif
!            get wz0=wz
            call cscalerProd(Tz, aTrack%w, w%z)
            temp = leng/rgt
            sinl = sin(temp)         ! this may be < 0
            cosl = cos(temp)

            w%x = sint*sinl
            w%y = sint*cosl
!///////////
!            call checkuv('w ', sint, w)
!///////////
            r%x = - rgt*sint*(cosl - 1.d0)
            r%y = rgt*sint*sinl
            r%z = leng*w%z
!                convert to original system.
            dr%x =
     *          Tx%x*r%x + Ty%x*r%y + Tz%x*r%z
            dr%y =
     *          Tx%y*r%x + Ty%y*r%y + Tz%y*r%z
            dr%z =
     *          Tx%z*r%x + Ty%z*r%y + Tz%z*r%z
!                directon cos
            dir%x =
     *          Tx%x*w%x + Ty%x*w%y + Tz%x*w%z
            dir%y =
     *          Tx%y*w%x + Ty%y*w%y + Tz%y*w%z
            dir%z =
     *          Tx%z*w%x + Ty%z*w%y + Tz%z*w%z
         else
            dr%x = 0.
            dr%y = 0.
            dr%z = 0.
            dir = aTrack%w
         endif
      else
         dr%x = leng*aTrack%w%x
         dr%y = leng*aTrack%w%y
         dr%z = leng*aTrack%w%z
         dir = aTrack%w
      endif
      end
