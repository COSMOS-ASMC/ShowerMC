!      Compute magnetic deflection (angle and displacement)
!      
!     ************************************************  
      subroutine cmagneticDef(aTrack, B, leng, dr,  dir)
!     ************************************************  
!    This is exact if B is const.
!   Note that:
!      dr must be added  as  newpos = oldpos +  dr
!      dir is the new direction cos.  These two are differenct
!      from the old version.
!
 
      implicit none
#include  "Ztrack.h"
! #include "Zmagfield.h"

      type(track)::aTrack       !input. a charged ptcl
      type(magfield)::B         !input. magnetic field vector.

!                         ptcl and magfiled coord is in Exyz

      real*8  leng                !input. length  travelled in m
      type(coord)::dr           !output. displacement vector in m
      type(coord)::dir          !output. new direction cos. 
!
      type(coord)::Tx
      type(coord)::Ty
      type(coord)::Tz
      type(coord)::w
      type(coord)::r

      real*8 Babs, sint, cost, rgt, sinl, cosl
      real*8   pabs, temp
!
!       get abs(B)
      Babs = sqrt( B%x**2 + B%y**2 + B%z**2 )
      if(Babs .gt. 0.) then
!            p in GeV
         pabs = sqrt(aTrack%p%fm%p(4)**2 - aTrack%p%mass**2)
         if(pabs .gt. 0.) then
!              transverse gyroradius in m. B in T, p in GeV
!              For negative charge, rgt < 0. This is o.k
            rgt = 3.3358*pabs/aTrack%p%charge/Babs

            Tz%r(1) = B%x/Babs
            Tz%r(2) = B%y/Babs
            Tz%r(3) = B%z/Babs

            call cscalerProd(aTrack%vec%w, B, cost)
            cost =max( min(cost/Babs, 1.0d0), -1.d0)
            sint = sqrt(1.0d0-cost**2)
            if(abs(sint) .gt. 1.d-6) then            
               call cvecProd(aTrack%vec%w, Tz, Tx)
               Tx%r(1) = Tx%r(1)/sint
               Tx%r(2) = Tx%r(2)/sint
               Tx%r(3) = Tx%r(3)/sint
               call cvecProd(Tz, Tx, Ty)
!////////////////
!               call checkuv('Tx', sint, Tx)
!               call checkuv('Ty', sint, Ty)
!               call checkuv('Tz', sint, Tz)
!               call checkorth('TxTy', sint, Tx, Ty)
!               call checkorth('TxTz', sint, Tx, Tz)
!//////////////////////               
            else
               Tx%r(1) = 1.d0
               Tx%r(2) = 0.
               Tx%r(3) = 0.
               Ty%r(1) = 0.
               Ty%r(2) = 1.d0
               Ty%r(3) = 0.
            endif
!            get wz0=wz
            call cscalerProd(Tz, aTrack%vec%w, w%r(3))
            temp = leng /rgt
            sinl = sin(temp)         ! this may be < 0
            cosl = cos(temp)

            w%r(1) = sint*sinl
            w%r(2) = sint*cosl
!///////////
!            call checkuv('w ', sint, w)
!///////////
            r%r(1) = - rgt*sint*(cosl - 1.d0)


            r%r(2) = rgt*sint*sinl
            r%r(3) = leng*w%r(3)
            
!                convert to original system.
            dr%r(1) =
     *          Tx%r(1)*r%r(1) + Ty%r(1)*r%r(2) + Tz%r(1)*r%r(3)
            dr%r(2) =
     *          Tx%r(2)*r%r(1) + Ty%r(2)*r%r(2) + Tz%r(2)*r%r(3)
            dr%r(3) =
     *          Tx%r(3)*r%r(1) + Ty%r(3)*r%r(2) + Tz%r(3)*r%r(3)


            dir%r(1) =
     *          Tx%r(1)*w%r(1) + Ty%r(1)*w%r(2) + Tz%r(1)*w%r(3)
            dir%r(2) =
     *          Tx%r(2)*w%r(1) + Ty%r(2)*w%r(2) + Tz%r(2)*w%r(3)
            dir%r(3) =
     *          Tx%r(3)*w%r(1) + Ty%r(3)*w%r(2) + Tz%r(3)*w%r(3)
         else
            dr%r(1) = 0.
            dr%r(2) = 0.
            dr%r(3) = 0.
            dir = aTrack%vec%w
         endif
      else
         dr%r(1) = leng*aTrack%vec%w%r(1)
         dr%r(2) = leng*aTrack%vec%w%r(2)
         dr%r(3) = leng*aTrack%vec%w%r(3)
         dir = aTrack%vec%w
      endif
      end
!c///////////////
!      subroutine checkuv(com, sint, v)
!      implicit none
!#include "Zcoord.h"
!      character*(*) com
!      record /coord/ v
!
!      real*8 eps, sint
!      eps = abs( (v.r(1)**2 + v.r(2)**2 + v.r(3)**2)-1.d0 )
!      if( eps .gt. 1.d-5) then
!         write(0,*) 'not unit', com, eps, sint
!      endif
!      end
!      subroutine checkorth(com, sint, v, u)
!      implicit none
!#include "Zcoord.h"
!      character*(*) com
!      record /coord/ v, u
!      real*8 temp, sint
!      call cscalerProd(u, v, temp)
!      if(abs(temp) .gt. 1.d-5) then
!         write(0,*) ' uxv nrt orth', com,  temp, sint
!      endif
!      end




