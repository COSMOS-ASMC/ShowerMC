!   If you want to modify the default treatment,  give AutoEmin in your
!   epicsfile a value >  1 and implement  a code below corresponding 
!   to that value. (if  0, no AutoEmin is used).

      subroutine epAutoEmin( gramage, media,
     *   EminG, EminE, Recoilmin )
!      This will be used when AutoEmin >= 1 in "epicisfile"
!      When a detector configuation is read, the minimum
!      energies for gamma, electron and recoil electron to be
!      sampled by random process, are fixed depending on the 
!      minimum thickness  of each component. 
!      Other minimum energies
!       charged particle (non-e) minimum.  made to be the same as gamma.
!       neutron kinetic energy.  fixed by epicsfile (5 MeV is default) 
!
!      The default treatment below may be modified by the user
!      if the default  below is not good for your applications.
!       
! **** Note: the treatment here is overriden by data in the "Modify" file
!      if the use of "Modify" file  is specified for that component in the
!      config file. 
      implicit none
#include "Zmass.h"
#include "Zmedia.h"
#include "ZepTrackp.h"

      real(8),intent(in):: gramage  ! "minimum" thickness of the component
                      !       in unit of g/cm2. (minimum is not necessary
                      !       minumum. It difficult to say it for, eg
                      !       prism.
       type(epmedia)::  media  ! input.  media of the component.
                      ! some vales for interest may be
                      ! media.X0  : r.l in cm
                      ! media.X0g :  //    in g/cm2
                      ! meida.I   : average ionization potential. in GeV.
      real(4),intent(out):: EminG, EminE, Recoilmin 
                      ! GeV). g, e, recoil min.  energy. 
                      ! Note.  EminE is total energy.
!***   For thin Si, e.g 32 micron m to 1.04 mm, all min=10keV is good
!***   for 1.04 micorn, 15keV is almost no change from 10 keV
!***   setting 3 keV for 32 micron  case is NG.

      !------------------ default ------------
      if( AutoEmin == 1 ) then  ! v9.154 use same one as AutoEmin =1, 2
!            set min between 10 keV to 100 keV.
         EminG = min( max(sqrt(gramage)*30.e-6, 10.e-6), 100.e-6)
         EminE = EminG + masele
         Recoilmin = max(EminG, 10.0*media%Z2*1.d-9)  ! For thin Pb,W
                                                 ! maybe some problems
      elseif( AutoEmin <= 3 ) then   ! 
         EminG = min( max(sqrt(gramage)*30.e-6, 10.e-6), 100.e-6)
         EminE = EminG + masele
         Recoilmin = max(EminG,  10.0*media%Z2*1.d-9)
    !--------------------------
      else
         write(0,*) ' AutoEmin=',AutoEmin, 
     * ' not yet suported; see prog/UserMayChange/epAutoEmin.f'
         stop
      endif
!      if( AutoEmin == 1 ) then ! old one; for reference kept
!!            set min between 10 keV to 100 keV.
!         EminG = min( max(sqrt(gramage)*100.e-6, 10.e-6), 100.e-6)
!         EminE = EminG + masele
!         Recoilmin = max(EminG, 14.0*media.Z2*1.d-9)
!
!      elseif( AutoEmin == 2 ) then  ! old 
!         EminG = min( max(sqrt(gramage)*150.e-6, 10.e-6), 150.e-6)
!         EminE = EminG*2 + masele
!         Recoilmin = max(EminG, 14.0*media.Z2*1.d-9)
!
      end subroutine epAutoEmin
