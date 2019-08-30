!  This version is to use older dpmjet3 implemented in Cosmos/Eics
!  result is same as the one in the upper diectoru (getXsec.inp)
#include "ZcosmosBD.h"
!$ CREATE DPMJET.FOR
!COPY DPMJET
!
!===program dpmjet=====================================================*
!
      PROGRAM DPMJET

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

! block data in DPMJET library (uncomment these declarations if library
! option is used)
!     EXTERNAL DT_BDEVAP,DT_BDNOPT,DT_BDPREE,DT_HADPRP,DT_BLKD46,
!    &         DT_BLKD47,DT_RUNTT,DT_NONAME,DT_ZK,DT_BLKD43

!     EXTERNAL PYDATA

! event flag
      COMMON /DTEVNO/ NEVENT,ICASCA

      character(120):: path

!-----------------------------------------------------------------------
! initialization

!   the following statement provides a call to DT_USRHIS(MODE=1) for
!   histogram initialization etc.
!      CALL DT_DTUINI(NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
      call cintModels('epics')
      call cformFullPath('CC.inp', path)
      call cinidpmjet(path)
      write(0,*) 'after--------------------'

!-----------------------------------------------------------------------
! generation of events

      DO 1 IEVT=1,NEVTS

!   some defaults, do not change!
         NEVENT = IEVT
         KKMAT  = -1
         ELAB   = EPN
!   uncomment if dpmjet3 is linked to particle transport code
!        ICASCA = 1

!***********************************************************************
! The following lines show how to select the target nucleus for runs
! with composite targets (and fixed projectile and energy!).
!
!   Sampling of the target nucleus (mass number NTMASS, charge NTCHAR)
!   according to the fractions defined with EMULSION input-cards.
!   The different nuclei are numbered as KKMAT = 1,2,3,...  according to
!   their appearance in the input-file.
         IF (IEMU.GT.0) THEN
!   Replace this selection by your own one if needed.
            CALL DT_GETEMU(NTMASS,NTCHAR,KKMAT,0)
!   Kkmat has to be negative for composite targets!
            KKMAT = -KKMAT
         ENDIF
!***********************************************************************

!***********************************************************************
! The following lines show how to define projectile, target and energy
! for this event in runs with Glauber-data file pre-initialized for a
! certain range of projectiles, targets and energies. The definitions
! have to be within the pre-initialized parameter range.
!
!   projectile-id (for hadron projectiles)
!        IDP    = 1
!   projectile mass and charge numbers
!        NPMASS = 12
!        NPCHAR = 6
!   target mass and charge numbers
!        NTMASS = 16
!        NTCHAR = 8
!   lab energy
!        ELAB = 200.0D0
!***********************************************************************

!***********************************************************************
! If an energy-range has been defined with the ENERGY input-card the
! laboratory energy ELAB can be set to any value within that range. For
! example:
!        ELO  = 10.0D0
!        EHI  = 1000.0D0
!        ELAB = DT_RNDM(ELAB)*(EHI-ELO)+ELO
!***********************************************************************

!   sampling of one event
         CALL DT_KKINC(NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,ELAB,KKMAT,IREJ)
         IF (IREJ.NE.0) GOTO 1

!   the following statement provides a call to DT_USRHIS(MODE=2) from
!   where the final state particles can be obtained

         CALL PHO_PHIST(2000,DUM)

    1 CONTINUE

!-----------------------------------------------------------------------
! output, statistics etc.

!   the following statement provides a call to DT_USRHIS(MODE=3) in
!   order to calculate histograms etc.
      CALL DT_DTUOUT

      END
