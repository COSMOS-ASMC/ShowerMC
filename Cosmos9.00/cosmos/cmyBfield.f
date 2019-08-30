!     This is the template for a non-Earth magnetic field.
!     If you want to give ObjFile /= " " to specify non-Earth object,
!     this file must be copied to the user's app. folder,
!     and modify the program to give a magnetic field (MagF)  at a given
!     location (pos).
!**** The user must comment out the next line at the top part of  this sub.****
!         call ccheckBfieldSub
!     If the user want to use a similar formula as geomagnetic one, 
!     Tracking/Geomag/cgeomag.f  may be copied and modified to make a new mag
!     field.   
!      
      subroutine cmyBfield( yearin, pos, MagF, icon )
      use modAtmosDef
      implicit none
#include "Ztrackp.h"
#include  "Zglobalc.h"
#include  "Zcoord.h"
#include  "Zmagfield.h"

      real(8),intent(in):: yearin ! e.g 2020.5; for non-Earth, probably not used
      type(coord),intent(in)::pos !  position data where MagF is to be computed
!      If we follow the case of cgeomag, pos and and MagF would be as follows:
!       If pos%sys is not "llh", pos is converted to "llh" system and
!       saved in  "cdata" (see below), and  it is used.  
      type(magfield),intent(out)::MagF ! magnetic field vector in ned;
!     ned means  MAgF =  (Bnorth, Beast, Bdown).  Unit is T.
!     It could be given in "xyz" syste with MagF%sys = "xyz"
!     When used in the tracking, 'xyz' system value is used.

!     The input/output arguments are the same as standard geomagnetic
!     calculation routine: cgeomag, though some of them may not be needed
!     for this case.      
!   
!        To convert coordinate system, you may use
!     call ctransCoord2( sys, inpos, outpos)
!        sys is target coord sys to be used in outpos
!     one of  "xyz", "llh", "sph"
!     xyz:  Earth center system; E-xyz. X is directed to  longitude 0, latitude 0
!                  Y is  90 deg. east, Z to the north. (m)
!     llh:  (latitude(deg), longitude(deg), height(m)) .
!     sph:  Polar coord sys. (teta(deg), fai(deg), r(m))
!     Another routine
!        call  ctransMagTo( sys, pos, a,  b)
!     sys: target coord. sys. of b  ( one of  'xyz', 'hva', 'ned')
!     pos: position where magnetic field a is given.
!     a: input B field. a%sys must be one of  'xyz', 'hva', 'ned'
!     b: output B field.   b%sys = sys
!      
!     inpos is the input pos; outpos is the ouput pos.

      integer,intent(out):: icon ! 0 if  ok


!
!
!       real*4 r, sumn, sume, sumd, t, cost, sint, x, tlonr, gmnc, 
!     *       cosml, sinml, hmnc, temp
!       real*8 gn, ge, gd
!       real*4 ssumd, ssumn, ssume 
!     integer m, n
      
       type(coord)::cdata

! this line must be comment out for actual application
       call ccheckBfieldSub
       
!
!     If "llh" system need not be used and the coord system is always fixed,
!        next few lines could be dropped or modified.
       if( pos%sys .eq. 'llh') then
          cdata = pos
       else     ! convet to llh
          call ctransCoord2('llh', pos, cdata)
       endif   
!
       icon =  0
!        if Bx, By, Bz in NED system are obtained 
!     call csetMagField('ned', Bx,By,Bz, MagF). If 'xyz' system,
!     call csetMagField('xyz', Bx,By,Bz, MagF). If 'xyz' system,

       end
      subroutine ccheckBfieldSub
      implicit none
      write(0,*)
     * 'You have given a non blank data to "ObjFile" to spcify'
      write(0,*)
     * '  non-Erath environment but '
      write(0,*)      
     * ' a) You might forget to copy $COSMOSTOP/cosmos/cmyBfield.f '
      write(0,*)
     * '  to your application folder and/or compile it together with '
      write(0,*)
     * '  your applicaion  OR'
      write(0,*)
     * ' b) forget to comment out the line "call ccheckBfieldSub"'
      write(0,*)
     *     ' in the top part of "cmyBfield" '
      stop
      end
      
