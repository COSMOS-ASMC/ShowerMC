      subroutine cmyEfield(aTrack, Efout)
!      This is a template for giving electric field which is
!      complex and   cannot  be expressed by the default method.
!      This should be copied to your application area and
!     modified as you like; in such a situation, the copied and
!     modified one  will override this one by usual "make".
!   
!     For actual run, you must specify HowEfield=2 in the
!     namelist ($HPARAM part).
!     If you don't prepare a copied one, this template is
!     compliled to avoid "missing subroutine" but if
!     HowEfield=2 is spcicified in such a case, stop will
!     be made:
!     Therefore, you must comment out the first exectutable
!     in this subroutine,
!           call ccheckEfieldSub      
!
      implicit none
#include "Zglobalc.h"
#include "Ztrack.h"
#include "Zobs.h"
#include "Zobsv.h"
#include "Zincidentv.h"
      type(track):: aTrack  ! input. Track info in E-xyz.  
       ! position and "where" info (or others ) are used.      
      real(8),intent(out):: Efout(3) ! Elecric field vector in
                 ! the E-xyz system. unit is  V/m.
!   This following information may be used to fix      Efout
!
!   particle position.
!    aTrack.pos.height: vertical height in m from the sealevel
!    aTrack.pos.depth:  depth in kg/m2.
!    aTrack.pos.xyz.r(1:3):  xyz pos in E-xyz system.
!                   to convert to other system,
!         see cdet2xyzFE.f and cxyz2det.f in Cosmos/Tracking.
!         Use routines cxyz2***.  *** depends on the requirment.
!         
!    aTrack.where: the observation level to which the ptcl
!                  is directed.     
!  The next observation point
!       ObsSites(aTrack.where).pos.height  ! height in m
!       ObsSites(aTrack.where).pos.depth  ! height in m
!  Deepest observation point
!       ObsSites(NoOfSites).pos.height
!       ObsSites(NoOfSites).pos.depth
!  distance r to the shower axis (perpendicular to the axis). r in m.
!       real(8)::r
!       call cgetPcoreDist(aTrack, r) 
!  distance r to the shower axis (horizontal to the axis). r in m
!       call cgetHcoreDist(aTrack, r)
!  time information 
!       aTrack.t/c*Tonsec  in ns. 
!      (c and Tonsec is deifined in Zglobalc.h; 0 is first col. 
!         point)
!  to convert Eelectric field (Ef) in detector or primary system
!      to the Exyz system. 
!     real(8):: Ef(3) 
!     call cdet2xyzD(Ef, Efout)  ! for detector to Exyz
!     call cprim2xyzD(Ef,Efout)  ! for primary to Exyz  
      real(8)::h
      real(8)::x,y
      real(8):: Efz, Ef(3)
      type(coord):: xyz
!     ****     comment out next when you modify the copied one ****
      call ccheckEfieldSub

      h = aTrack.pos.height
      
      if(h > 7e3 ) then
         Efz = 0.
      else
!           convert x,y. origin is  at the deepest det. pos
         call cxyz2det( ObsSites(NoOfSites).pos.xyz,
     *            aTrack.pos.xyz, xyz)
         x = xyz.r(1)
         y = xyz.r(2)
!             E eists only in  some area
         if( x > 0. .and. y > 0.) then
            Efz= 360d3*sin(3.1415/2.*h/2d3)*exp(-h/5.d3)
         else
            Efz = 0.
         endif
      endif
      Ef(:) =(/0.d0, 0.d0, Efz/)
      call  cdet2xyzD(Ef, Efout)
      end subroutine cmyEfield

      subroutine ccheckEfieldSub
      implicit none

      write(0,*)
     *  'Did you copy $COSMOSTOP/cosmos/cmyEfield.f to your app. '
      write(0,*)
     *  ' folder and modified it ?'
      write(0,*)
     * ' If so, fogot to comment out a line: "call ccheckEfieldSub" ? '
      write(0,*)
     * '    or forgot "make clean ;make" after making modified one ?'

      
      write(0,*)       
     * ' Or You have given "HowEfield=2" but not yet made '
      write(0,*)
     *     ' a copy of cmyEfield ? '

      stop
      end
      
