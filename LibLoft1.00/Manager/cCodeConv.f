      module modCodeConv
      implicit none

      interface  cpdg2cos
      module procedure cpdg2cosPtcl
      module procedure cpdg2cosI
      module procedure cpdg2cosHI      
      end interface

      interface  ccos2pdg
      module procedure ccosPtcl2pdg
      module procedure ccosI2pdg
      module procedure ccosHI2pdg
      end interface

      contains


      subroutine cpdg2cosPtcl(pdgc, aPtcl)
      implicit none
#include "Zptcl.h"
      integer,intent(in):: pdgc
      type(ptcl),intent(out):: aPtcl

      call ckf2cosB(pdgc, aPtcl)

      end       subroutine cpdg2cosPtcl

      subroutine cpdg2cosI(pdgc, code, subcode, charge)
      implicit none
      integer,intent(in):: pdgc
      integer,intent(out):: code, subcode, charge
      call ckf2cos(pdgc, code, subcode, charge)
      end       subroutine cpdg2cosI

      subroutine cpdg2cosHI(pdgc, code, subcode, charge)
      implicit none
      integer,intent(in):: pdgc
      integer(2),intent(out):: code, subcode, charge
      integer:: cd, sc, cg
      call ckf2cos(pdgc, cd, sc, cg)
      code= cd
      subcode = sc
      charge = cg
      end       subroutine cpdg2cosHI
!   **************
      subroutine ccosPtcl2pdg(aPtcl, pdgc)
      implicit none
#include "Zptcl.h"
      type(ptcl),intent(in):: aPtcl
      integer,intent(out):: pdgc


      call ccos2kfB( aPtcl, pdgc)
      end  subroutine ccosPtcl2pdg

      subroutine ccosI2pdg(code, subcode, charge, pdgc)
      implicit none
      integer,intent(in):: code, subcode, charge
      integer,intent(out):: pdgc
      
      call ccos2kf( code, subcode, charge, pdgc)
      end subroutine ccosI2pdg

      subroutine ccosHI2pdg(code, subcode, charge, pdgc)
      implicit none
      integer(2),intent(in):: code, subcode, charge
      integer,intent(out):: pdgc

      integer:: cd, sc, cg
      
      cd = code
      sc = subcode
      cg = charge
      call ccos2kf( cd, sc, cg, pdgc)
      end subroutine ccosHI2pdg
      
      end module modCodeConv
