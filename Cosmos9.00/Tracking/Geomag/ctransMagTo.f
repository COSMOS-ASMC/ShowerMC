!
!       transform magnetic field components in one coordinate
!       sytem to another.
!  a.sys     \sys      'xyz'     'hva'   'ned'
!  
!      'xyz'             o         o        o
!      'hva'             o         o        o
!      'ned'             o         o        o 
!
      subroutine ctransMagTo(sys, pos, a, b)
!
!      sys: character*(*).  input.  'xyz', 'hva', 'ned'
!                           the target coordinate system where
!                           magnetic filed is represented.
!      pos: /coord/.        input.  position where mag is given
!        a: /magfield/  input. 
!        b: /magfield/  output. transformed component, b.sys=sys
!     
      implicit none

#include  "Zcoord.h"
#include  "Zmagfield.h"
      character*(*) sys
      type(magfield)::a
      type(magfield)::b
      type(coord)::pos
!
      character*70  msg
!
      if( a%sys .eq. sys) then
         b = a
      elseif(sys .eq. 'xyz') then
         call cMag2eCent(pos, a, b) 
      elseif(a%sys .eq. 'xyz') then      ! current system
         if(sys .eq. 'hva') then
             call cxyz2hva(pos, a, b)
         elseif(sys .eq. 'ned') then
             call cxyz2ned(pos, a, b)
         else
             write(msg, *) ' ctransMagTo: sys=', sys, ' invalid'
             call cerrorMsg(msg, 0)
         endif
      elseif(a%sys .eq. 'hva') then
         if( sys .eq. 'ned') then
                call chva2ned(a, b)
         else
             write(msg, *) ' ctransMagTo: sys=', sys, ' invalid'
             call cerrorMsg(msg, 0)
         endif
      elseif(a%sys .eq. 'ned') then
         if(sys .eq. 'hva') then
             call cned2hva(a, b)
         else          
             write(msg, *) ' ctransMagTo: sys=', sys, ' invalid'
             call cerrorMsg(msg, 0)
         endif
      else
         write(msg, *) ' ctransMagTo: a%sys=', a%sys, ' invalid'
         call cerrorMsg(msg, 0)
      endif
      end
!------------------------------------------------------------- 
      subroutine cMag2eCent(pos, a, b)
!               to earth_center system (xyz system)
        implicit none
#include  "Zcoord.h"
#include  "Zmagfield.h"
        type(magfield)::a
        type(magfield)::b
        type(coord)::pos
!
        type(coord)::postemp
        character*70  msg
!
        call ctransCoord2('llh', pos, postemp)
        if(a%sys .eq. 'ned') then
               call cned2eCent(postemp, a, b)
        elseif(a%sys .eq. 'hva') then
               call chva2ned(a, b)
               call cned2eCent(postemp, b, b)
        else
               write(msg, *) 'cMag2eCent: mag system=', 
     *         a%sys, 
     *         ' not yet supported'
               call cerrorMsg(msg, 0)
        endif   
       end
!----------------------------------------------------------------------
       subroutine cned2eCent(pos, a, b)
!          pos: /coord/  input.  pos.sys should be 'llh'
!            a: /magfield/ input.  in 'ned' system
!            b: /magfield/ output. in 'xyz'  system
!                b can be the same one as a.
       implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
!
       type(coord)::pos
       type(magfield)::a
       type(magfield)::b
       character*70 msg
!
       real*8 cosphi, sinphi, coslam, sinlam, x, y, z
!
       if(pos%sys .ne. 'llh') then
          write(msg, *)'cned2eCent: input pos%sys=',pos%sys,
     *              ' invalid. should be llh'
          call cerrorMsg(msg, 0)
       endif
#ifdef UNIONMAP
             cosphi = cos(pos%lat*Torad)
             sinphi = sin(pos%lat*Torad)
             coslam = cos(pos%long*Torad)
             sinlam = sin(pos%long*Torad)
             x = - (a%d *cosphi + a%n*sinphi) * 
     *           coslam 
     *                     - a%e *sinlam
             y = - (a%d *cosphi + a%n*sinphi)* sinlam
     *                     + a%e*coslam                  
             z = - a%d *sinphi + a%n* cosphi
#else
             cosphi = cos(pos%r(1)*Torad)
             sinphi = sin(pos%r(1)*Torad)
             coslam = cos(pos%r(2)*Torad)
             sinlam = sin(pos%r(2)*Torad)
             x = - (a%z *cosphi + a%x*sinphi) * 
     *           coslam 
     *                     - a%y *sinlam
             y = - (a%z *cosphi + a%x*sinphi)* sinlam
     *                     + a%y*coslam                  
             z = - a%z *sinphi + a%x* cosphi
#endif      
             call csetMagField('xyz', x, y, z, b)
        end
!------------------------------------------------------------
        subroutine cxyz2ned(pos, a, b)
        implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
!
        type(coord)::pos
        type(magfield)::a
        type(magfield)::b
        real*8 cosphi, sinphi, coslam, sinlam, x, y, z
        real*8  adcans, n, e, d
        character*70 msg
!
       if(pos%sys .ne. 'llh') then
          write(msg, *)'cxyz2ned: input pos%sys=',pos%sys,
     *              ' invalid. should be llh'
          call cerrorMsg(msg, 0)
       endif
       if(a%sys .ne. 'xyz') then
          write(msg, *) 'cxyz2ned: a%sys=', a%sys, ' invalid'
          call cerrorMsg(msg, 0)
       endif
#ifdef UNIONMAP
!
       cosphi = cos(pos%lat*Torad)
       sinphi = sin(pos%lat*Torad)
       coslam = cos(pos%long*Torad)
       sinlam = sin(pos%long*Torad)
       x = a%x
       y = a%y
       z = a%z
!        -(a.d*cosphi + a.n*sinphi) 
       adcans = x * coslam + y *sinlam
       d =-  ( adcans*cosphi + z * sinphi )
!      n = -x*sinphi*coslam - y*sinphi*sinlam + z*cosphi
       n = - adcans*sinphi + z * cosphi
!       e =  (y + (d*cosphi + n*sinphi)*sinlam) * coslam  -
!     *      (x + (d*cosphi + n*sinphi)* coslam) * sinlam
       e =  -x*sinlam + y*coslam
#else
!
       cosphi = cos(pos%r(1)*Torad)
       sinphi = sin(pos%r(1)*Torad)
       coslam = cos(pos%r(2)*Torad)
       sinlam = sin(pos%r(2)*Torad)
       x = a%x
       y = a%y
       z = a%z
!        -(a.z*cosphi + a.x*sinphi) 
       adcans = x * coslam + y *sinlam
       d =-  ( adcans*cosphi + z * sinphi )
       n = - adcans*sinphi + z * cosphi
!       e =  (y + (d*cosphi + n*sinphi)*sinlam) * coslam  -
!     *      (x + (d*cosphi + n*sinphi)* coslam) * sinlam
       e =  -x*sinlam + y*coslam
#endif      
       call csetMagField('ned', n, e, d, b)
       end
!------------------------------------------------------------
        subroutine cxyz2hva(pos, a, b)
        implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
!
        type(coord)::pos
        type(magfield)::a
        type(magfield)::b
!
        call cxyz2ned(pos, a, b)
        call cned2hva(a, b)
        end
      subroutine cned2hva(a, b)
!        transform magnetic components from norht-east-down system
!        to horizontal-vertical-deflection_angle system.
!      a:  /magfield/ input.
!      b:  /magfield/ output. b can be the same entity as a.
!             b.h: horizontal component
!             b.v: vertical component
!             b.a: deflection angle (deg).  + is from the north to
!                   the clockwise direction.
      implicit none
#include  "Zglobalc.h"
#include  "Zmagfield.h"
!
      type(magfield)::a
      type(magfield)::b
      real*8 h, v, ang
!
#ifdef UNIONMAP
      h = sqrt(a%n**2+a%e**2)
      if(a%e .eq. 0. .and. a%n .eq. 0.) then
         ang = 0.
      else
         ang = atan2(a%e, a%n)*Todeg
      endif
      v = a%d
#else
      h = sqrt(a%x**2+a%y**2)
      if(a%y .eq. 0. .and. a%x .eq. 0.) then
         ang = 0.
      else
         ang = atan2(a%y, a%x)*Todeg
      endif
      v = a%z
#endif
!
      call csetMagField('hva', h, v, ang, b)
      end
!---------------------------------
      subroutine chva2ned(a, b)
!          inverse of the above
      implicit none
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"
      type(magfield)::a
      type(magfield)::b
      real*8 n, e, d
!
#ifdef UNIONMAP
      n = a%h *cos(a%a*Torad)
      e = a%h *sin(a%a*Torad)
      d = a%v
#else
      n = a%x *cos(a%z*Torad)
      e = a%x *sin(a%z*Torad)
      d = a%y
#endif
      call csetMagField('ned', n, e, d, b)
      end
