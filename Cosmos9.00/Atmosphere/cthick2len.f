!
!     Thickness of air is converted into length.
!     This should be placed rather in Tracking directory.
!     ****************************************
      subroutine  cthick2len(aTrack, tin, leng, t, jcut)
!     ****************************************
!
!    aTrack;  considered track
!   tin: input. real*8.  thickness of air to be travelled in kg/m2.
!  leng: output. real*8. if jcut=0; length in m corresponding to tin (= t). 
!       if jcut=1; lenght in m correslponding to t ( < tin)
! t: output. real*8. if jcut=0, this becomes the same as tin
!                       else  a shorter value than tin  (kg/m2).  This is
!                       the thickness corresponding to leng (m).  
!               for neutrino case, this is not true.
! jcut: output. integer.  if 0, thickness tin is successfully converted into leng.
!                         else, tin is too long. It is cut to a value 't'
!                         so that the resultant accuracy of leng is high.
!                         cut is made if tin is longer than the (current vertical depth)/5.
!
!
      use modAtmosDef
      implicit none
!----      include 'Zearth.h'
! #include  "Zearth.h"
! #include  "Zatmos.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zcode.h"
      type(track)::aTrack

      real*8 z, cosz, tin, leng, t
      integer jcut
!
!    z: input. real*8. current vertical height in m
!  cosz:input. real*8. cos of the current zenith angle

      real*8 maxkgm2/100./

      real*8 cvh2den, s, cnewcos, cnewh, zp, cosp
      real*8 cutf/5./, cthick2h
!       Method.
!
!       If we expess a heigh z in terms of the distance from the
!    center of the earth, r = R + z (R: radius), and, the zenith angle
!    cos, 
!    r cos - r'cos' = s  (s is the distance between two points)
!    r sin = r'sin'
!    r^2 + s^2 -2rscos = r'^2
!   Hence,
!   z' = sqrt(r^2 + s^2 - 2rscos) - R
!      = sqrt( (r-scos)^2 +s^2(1-cos^2) ) - R
!      = h -s cos s^2(1-cos^2)/(r - scos)/2 + O(s^4)      
!      = h - s cos + s^2(1-cos^2)/2/r + s^3/r^2/2 cos(1-cos^2) + O(s^4)
!  The slant thickness in the length s is
!   
!     t = int[0,s](rho(s)ds
!       = int[0,s](rho(z + f(s, z))ds
!  where
!      f(s, z) = -s cos + s^2(1-cos^2)/r/2 + s^3/r^2/2 cos(1-cos^2)
!
!    rho(z + f(s,z)) = rho(z) + rho'(z)f(s,z) + rho''(z)/2 f(s,z)^2
!
!  Taking upto s^3, the indefinit integral is
!
!   f1(s,z) =int(f(s,z)ds = -s^2/2 cos + s^3/r/6(1-cos^2) 
!
!   f2(s,z) =int(f^2)ds = s^3/3 cos^2
!
! 
!
      real*8 f1, f2, int1, int2,  rp  !  , clenbetween2h
      real*8 cvh2denp, cvh2den2p, rho, rhop, d, cvh2thick
      real*8 ct2lTA

      real*8 nearv/0.5/, nearv2/0.99/, t1, t2
      real*8 maxact
      integer i, icon, iterate

      f1(s) = -s**2/2*cosp + s**3/rp/6.d0 * (1.d0-cosp**2)
!     *    +  s**4/rp/rp/8.d0*cosp*(1.d0-cosp**2)    ! negligible

      f2(s) = s**3/3.d0 * cosp**2
!     *     - s**4/4.d0/rp*cosp*(1.d0-cosp**2)   ! negligible

!----------------------------------------------------------------------
      z = aTrack%pos%height
      cosz = aTrack%vec%coszenith
      if( aTrack%p%code .eq. kneumu .or.
     *    aTrack%p%code .eq. kneue ) then
         leng = 1.e5
         jcut = 1
         t = tin
         return   !   **************
      endif
      if( tin .eq. 0.) then
         leng=0.
         jcut = 0
         t = tin
         return  ! *********
      endif

      if(UseTbl .and. z .lt. Htop ) then
         maxact = 100.d3
      elseif(abs(cosz) .gt. nearv) then
         maxact = 20*maxkgm2 *(abs(cosz) - nearv) + maxkgm2
      else
         maxact = maxkgm2
      endif
      t = min(tin, maxact)
      if(t .eq. maxact) then
         jcut= 1
      else
         jcut = 0
      endif

      if(.not. (UseTbl .and. z .lt. Htop) .and. 
     *    (abs(cosz) .lt. nearv) ) then
!               max movable thickness
         t = min(t, cvh2thick(z)/cutf/(abs(cosz) + .1) )
         if(t .ne. tin) then
           jcut =1
        endif
      endif
!     
      iterate =0   ! flag for later use. 

      if(UseTbl .and. z .lt. Htop) then
         leng = ct2lTA(z,  t)
!      elseif( abs(cosz) .gt. nearv .or. 
!     *  ( z .lt. 15.d3  .and. t .lt. 350.0d0 ) ) then
      elseif( abs(cosz) .gt. nearv ) then
         t1 = cvh2thick(z)
         t2 = t1 + cosz* t
         if(t2 .le. 0.) then
            t2 = t1/10.
            t = (t2 - t1)/cosz
            jcut = 1
         endif
         zp = cthick2h(t2)
!         leng = clenbetween2h(z+Eradius, zp+Eradius, cosz)
         call clenbetw2h(z+Eradius, zp+Eradius, cosz, leng, icon)

         
         if(icon .ne. 0) then
!            cannot reach to zp. This may happen if z is very large
            iterate = 1
         else
!            if(abs(cosz) .lt. nearv2) then
            if( t .gt. 350.d0 ) then
!              correction by using  cos at middle of leng
               cosp = cnewcos(z+Eradius, cosz, leng/2)
               t2 = t1 + cosp * t
               if(t2 .le. 0.) then
!                   we shall move by real length rather than by thickness
                  jcut = 1
                  leng = leng/2
                  zp = cnewh(z+Eradius, cosz, leng) - Eradius
                  t2 = cvh2thick(zp)
                  t =( t2 - t1 )/((cosz+cosp)/2)
               else
                  zp = cthick2h(t2)
!                  leng = clenbetween2h(z+Eradius, zp+Eradius, cosz)
                  call clenbetw2h(z+Eradius, zp+Eradius, cosz,
     *                           leng, icon)
                  if(icon .ne. 0) then
                     iterate = 1
                  endif
               endif   
            endif
         endif
      else         
         iterate = 1
      endif
      
      if( iterate .ne. 0) then
         rho = cvh2den(z)
         rhop = cvh2denp(z)
!               at very high altitude,  rhop = 0
         if(cosz .le. 0.d0 .or. rhop .le. 1.0d-18)  then
            s = t/rho           !  length if density is const.
         else
            d = rho**2 - 2* rhop *cosz * t
            s = (rho -sqrt(d) )/ rhop/cosz
         endif

      
         if(rhop .gt. 1.0d-18) then
            do i=1,3
!               get height and cos at s/2 ahead.
               cosp = cnewcos(Eradius+z, cosz, s/2)
               zp = cnewh(Eradius+z, cosz, s/2) - Eradius
!                 write(0,*) " cosp =", cosp, " height at z-s/2*cos=",zp
!                 once again
!                s = t/cvh2den(zp)   !  length if density is const.
!                 cosp = cnewcos(Eradius+z, cosz, s/2)
!                zp = cnewh(Eradius+z, cosz, s/2) - Eradius
!
!               expand rho at zp, and get integrals of f and f^2
!                 [ -s/2 to s/2] ds
!
               rp = zp + Eradius
               int1 = f1(s/2) - f1(-s/2)
               int2 = f2(s/2) - f2(-s/2)
!                  get new s
               s = (t - cvh2denp(zp)*int1 - cvh2den2p(zp)*int2/2)
     *                 /  cvh2den(zp) 
            enddo
         endif
         leng = s
      endif
!  &&&&&&&&&
      if(leng .le. 0.) then
         rho = cvh2den(z)
         leng =t / rho
      endif
      end

      subroutine cvh2mediaNo(vh, mno)
      use modAtmosDef
      implicit none
!  #include "Zatmos.h"
#include "ZmediaLoft.h"
      real(8),intent(in):: vh   !  vertical height in m
      integer,intent(out):: mno ! media No

      integer:: loc

      if( NoOfMedia == 1 ) then
         mno = 1
      else
         call kdwhereis(vh, atmos%nodes, atmos%z, 1,  loc)
         if(loc == 0 ) then
            loc = 1
         elseif(loc > atmos%nodes) then
            loc = atmos%nodes
         endif
         mno = atmos%node2mediaNo(loc)
      endif
      end   subroutine cvh2mediaNo

      subroutine cvh1vh2Mcheck(vh1, vh2, node)
! see if media at vh1 to vh2 are the same. then, node =0
!     if media changes at some node, node > 0 is given.
      use modAtmosDef
      implicit none
#include "ZmediaLoft.h"
      
      real(8),intent(in):: vh1  ! vertial height  in m
      real(8),intent(in):: vh2  ! vertial height  in m   vh2 need not be   > vh1.
      integer,intent(out):: node
      
      integer::i, loc1, loc2

      if( NoOfMedia == 1 ) then
         node= 0
      else
         call kdwhereis(vh1, atmos%nodes, atmos%z, 1,  loc1)
         call kdwhereis(vh2, atmos%nodes, atmos%z, 1,  loc2)
         node = 0
         if(loc2 > loc1 ) then
            do i = loc1, loc2-1
               if( atmos%node2mediaNo(i)
     *              /=  atmos%node2mediaNo(i+1) ) then
                  node=i+1
                  exit
               endif
            enddo
         elseif( loc1 > loc2 ) then
            do i = loc1, loc2+1, -1
               if( atmos%node2mediaNo(i)
     *              /=  atmos%node2mediaNo(i-1) ) then
                  node = i-1
                  exit
               endif
            enddo
         endif
      endif
      end   subroutine cvh1vh2Mcheck

      subroutine cnearestDiffMedia(node, nodeB, nodeU)
      use modAtmosDef
      implicit none
#include "ZmediaLoft.h"      
      integer,intent(in):: node  ! reference node
      integer,intent(out):: nodeU !   the  nearest node # (> node) where
!     media is /= current media
!        If the diff. media does not exist,  0  is put
      integer,intent(out):: nodeB !   the  nearest node # (< node) where
!              media is /= current media
!        If the diff. media does not exist,  0 is put

      integer:: i, mno

      nodeB = 0
      nodeU = 0
      if( NoOfMedia > 1 )  then
         mno=atmos%node2mediaNo(node)
         do i = node+1, atmos%nodes
            if( atmos%node2mediaNo(i) /= mno) then
               nodeU=i
               exit
            endif
         enddo
         do i = node -1, 1, -1
            if( atmos%node2mediaNo(i) /= mno) then
               nodeB = i
               exit
            endif
         enddo
      endif
      end       subroutine cnearestDiffMedia
