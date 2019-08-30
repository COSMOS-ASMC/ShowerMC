!     ****************************************************************
!     *                                                              *
!     * epPhotoE:  to compute coefficient which appear in 
!     *            photoelectric effect.
!     *     must be called after epGetEffZA is called                *
!
!
!
      subroutine epSetPhotoE(media, pe)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"

       type(epmedia):: media  ! input. 
      intent(in):: media
       type(photoE):: pe    ! output. photo electric effect const
      intent(out):: pe
!
      real*8 az, z
      real*8 eg, prob, path, simple

      data eg/0.2d-3/
     
!      z = media.Zeff
      z = media%Z
      az = z/137.0

      if(z .lt. 31.0) then
         pe%b0 = 0.5/(  ( .9663*az + 5.023) *az + .9211)
      else
         pe%b0 = 1./(  ( .9663*az + 5.023) *az + .9211)
      endif
      pe%b1 = (2.56*az - 2.632)*az + 1.90
      pe%b2 = ( (-6.563*az+ 8.25)*az - 5.616)*az + 2.097
      pe%fa = ( .2762*az - .0288) *az + 1.083
      pe%a = masele - 13.6d-9*z**2     ! Kshell work func
      pe%b = masele - 3.4d-9*z**2      ! Lshell work func
      pe%p = 0.6 * media%Z5byAeff/137.0**4 *  media%X0g
      pe%ek = 13.6d-9*z**2 
      pe%l = (8. -z*2.5/80.)
      pe%cr = 1.0    ! tentative correction value
!////////////
!      write(0,*) ' in Set: media.pe.cr=',media.pe.cr
!///////////////
      if( z .lt. 31.0 ) then
!          get  cross section at eg = 0.2 MeV
!         call epphotoEp(pe, eg, prob, path)  < v8.0
         call epphotoEp(media, eg, prob, path) 
!//////////
!         write(0,*)  'after epphotoEp: media.ep.cr=',
!     *   media.pe.cr,' pe.cr=',pe.cr
!         write(0,*) ' prob =', prob, ' path =',path
!////////////
         prob= prob/media%mbtoPX0 !  mb/atom
!                    sqrt(32)
         simple = 5.65685d0/(eg/masele)**3.5 * az**4*z *665.0d0
!///////////
!         write(0,*) 'in Set eg =',eg, 'az=',az, ' z=',z
!/////////////////
         pe%cr = simple/prob
      endif
      end
