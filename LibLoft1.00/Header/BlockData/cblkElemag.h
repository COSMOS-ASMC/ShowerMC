!         block data cblkElemag


!#include "Zelemagp.h"
!
!         data 
!     *	RecoilKineMinE /0.2d-3/,     !  in modEMcontrol  100 keV
!     *  KnockOnRatio  /1.0d0/,       ! 
!     *  AnihiE /1.e3/,            ! Eposi < 1 TeV, anihilation considered
!     *  X0 /365.667/,                ! radiation length of air in kg/m2
!     *  Ecrit /81.e-3/,             ! critical energy of air in GeV
!     *  MaxComptonE /1./,       ! compton is considered below 1 GeV
!     *  MaxPhotoE /1.e-3/,            ! above this, PhotoElectric effect neg.
!     *  MinPhotoProdE /153.e-3/,     ! below 153 MeV, no gp --> hadrons
!     *	Es /19.3d-3/,                 ! scattering const. Note, not 21 MeV.

!     *  PhotoProd /.false./,        ! gp--> hadrons not considered.
!     *  Moliere /2/,          !  Moliere scattering
!     *  ALateCor /1/     ! Angle-lateral correlation for multiple scat.
!     *  Excom1 /1.0d-3/,                ! if Eg <=Excom1,  XCOM data is used for
!     *  Knockon /.true./,            ! knockon is considered. Obsolete	   
                                    ! compton/p.e/coherent .
!     *  Excom2 /1.0d0/,                ! if Eg <= Excom2, XCOM data for Pair cre.

 !        end
