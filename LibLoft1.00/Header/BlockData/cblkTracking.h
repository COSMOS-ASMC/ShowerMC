!     *************************
      block data  cblkTracking
!     *************************
      implicit none

#include  "Ztrackp.h"

      data 
     *  ExactThick /.false./ ,
     *  Freec /.true./ ,
     *  TraceDev /21/ , 
     *  TraceDir /' '/ ,
     *  BorderHeightH /0./ ,
     *  BorderHeightL /0./ ,
     *  HowGeomag /11/ ,
!     *  MagN /0./ ,
!     *  MagE /0./ ,
!     *  MagD /0./ ,
     *  RatioToE0 /1.e-5/ ,
     *  UseRungeKutta /0/ ,
     *  TimeStructure /.true./ ,
     *  Trace /0/ ,
     *  WaitRatio /0.d0/ ,	
     *  Truncc/5.0/ ,
     *  Truncn /1.e-2/ ,
     *  Truncx /2./
!     *  MagChgDist /20.e3/ ,
      data
!     * IncMuonPolari /.true./ ,
     * KEminObs /2*500.0d-6, 6*50.0d-3/ ,
     * ThinSampling /.false./,
     * EthinRatio /3.e-6, 10.0e5, 0., 0./ ,
     * Generate /'em'/ , 
     * Generate2 /' '/ ,
     * BackAngLimit  /-1.0/ ,
     * OneDim /0/
!     * MagBrem /2/,
!     * MagPair /1/,
!     * LpmEffect /.true./,
!     * MagPairEmin /2.e10/,         ! Eg > this ==> magnetic pair considered.> 3x10^19 eV
!     * MagBremEmin /3.e8/,          ! Ee > this ==> magnetic brems considered > 10^18 eV
!     * UpsilonMin  /3.e-3/,         ! if upsilon is < this, no magnetic brem
!     * LpmPairEmin /1.e9/,         ! Eg > LpmPairEmin && LpmEffect ==> LPM effect for g
!     * LpmBremEmin /1.e8/           ! Ee > LpmBremEmin && LpmEffect ==> LPM effect for e

       data
     * Reverse /0/,
     * PathLimit /13000d4/,         !  ~ 20 x Eradius
     * EndLevel/0/,
!     * MuNI /2/, 
!     * MuBr /2/,
!     * MuPr /2/,
     * LamorDiv /5./,
!     * Eabsorb /0,0/,
     * StepControl /2./,
!     * HowPhotoP /1/,
     * PhitsXs /0/,  JamXs /0/, AAXsec/0/, JamFragment/1/,
     * ASRforDPM /'m'/
      end

