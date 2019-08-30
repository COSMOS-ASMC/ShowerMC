!       *******************************************************
!       *  ciniSPrim: initialize csampPrimary
!       *******************************************************
!   Initializes a table for sampling the primary energy
!         and species.  For data to be prepared for the table, see
!         sample.d in Data/Primary directory.
!   NOTE:  Init for angle sampling is not done here. It is 
!          managed in the Tracking directory.
!
        subroutine ciniSPrim(fn)
!          fn:  character string.  input. file name where primary spectrum
!                                         data is given.  
!                                         See sample.d in Data/Primary
!
        implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
#include  "Zprimaryv.h"
        character*(*) fn
!
!     
        call csetPrimTbl      ! make PrimaryIdTbl
        call csetEUnitTbl     ! make ErgUnitTbl
        call ciniSPrim0(Prim, fn)
        end
!       *******************************************************
!       *  ciniSPrim0: initialize csampPrimary for a given rec
!       *******************************************************
!   Initializes a table for sampling the primary energy
!         and species.  For data to be prepared for the table, see
!         sample.d in Data/Primary directory.
!
        subroutine ciniSPrim0(prm, fn)
!          fn:  character string.  input. file name where primary spectrum
!                                         data is given.  
!                                         See sample.d in Data/Primary
!
        implicit none
#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"

        character*(*) fn
        type(primaries):: prm
        integer icon
!
!     
        call copenf(TempDev,  fn, icon)
        call cskipComment(TempDev, icon)
        call crdPrimData(prm )
        close(TempDev)
        call cprocPrimDt(prm )
        end
!
!     ******************************************     
      subroutine csetPrimTbl
         implicit none

#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
!
        PrimaryIdTbl(1)%symb = 'gamma'
        PrimaryIdTbl(1)%code = kphoton
        PrimaryIdTbl(1)%subcode = 0
        PrimaryIdTbl(1)%charge = 0

        PrimaryIdTbl(2)%symb = 'photon'
        PrimaryIdTbl(2)%code = kphoton
        PrimaryIdTbl(2)%subcode = 0
        PrimaryIdTbl(2)%charge = 0

        PrimaryIdTbl(3)%symb = 'e'
        PrimaryIdTbl(3)%code = kelec
        PrimaryIdTbl(3)%subcode = regptcl
        PrimaryIdTbl(3)%charge = -1

        PrimaryIdTbl(4)%symb = 'e-'
        PrimaryIdTbl(4)%code = kelec
        PrimaryIdTbl(4)%subcode = regptcl
        PrimaryIdTbl(4)%charge = -1

        PrimaryIdTbl(5)%symb = 'electron'
        PrimaryIdTbl(5)%code = kelec
        PrimaryIdTbl(5)%subcode = regptcl
        PrimaryIdTbl(5)%charge = -1

        PrimaryIdTbl(6)%symb = 'e+'
        PrimaryIdTbl(6)%code = kelec
        PrimaryIdTbl(6)%subcode = antip
        PrimaryIdTbl(6)%charge = 1

        PrimaryIdTbl(7)%symb = 'positron'
        PrimaryIdTbl(7)%code = kelec
        PrimaryIdTbl(7)%subcode = antip
        PrimaryIdTbl(7)%charge = 1

        PrimaryIdTbl(8)%symb = 'mu-'
        PrimaryIdTbl(8)%code = kmuon
        PrimaryIdTbl(8)%subcode = regptcl
        PrimaryIdTbl(8)%charge = -1

        PrimaryIdTbl(9)%symb = 'mu+'
        PrimaryIdTbl(9)%code = kmuon
        PrimaryIdTbl(9)%subcode = antip
        PrimaryIdTbl(9)%charge = 1

        PrimaryIdTbl(10)%symb = 'pi+'
        PrimaryIdTbl(10)%code = kpion
        PrimaryIdTbl(10)%subcode = regptcl
        PrimaryIdTbl(10)%charge = 1

        PrimaryIdTbl(11)%symb = 'pi-'
        PrimaryIdTbl(11)%code = kpion
        PrimaryIdTbl(11)%subcode = antip
        PrimaryIdTbl(11)%charge = -1

        PrimaryIdTbl(12)%symb = 'pi0'
        PrimaryIdTbl(12)%code = kpion
        PrimaryIdTbl(12)%subcode = 0
        PrimaryIdTbl(12)%charge = 0

        PrimaryIdTbl(13)%symb = 'k+'
        PrimaryIdTbl(13)%code = kkaon
        PrimaryIdTbl(13)%subcode = regptcl
        PrimaryIdTbl(13)%charge = 1

        PrimaryIdTbl(14)%symb = 'k-'
        PrimaryIdTbl(14)%code = kkaon
        PrimaryIdTbl(14)%subcode = antip
        PrimaryIdTbl(14)%charge = -1

        PrimaryIdTbl(15)%symb = 'k0l'
        PrimaryIdTbl(15)%code = kkaon
        PrimaryIdTbl(15)%subcode = k0l
        PrimaryIdTbl(15)%charge = 0

        PrimaryIdTbl(16)%symb = 'k0s'
        PrimaryIdTbl(16)%code = kkaon
        PrimaryIdTbl(16)%subcode = k0s
        PrimaryIdTbl(16)%charge = 0

        PrimaryIdTbl(17)%symb = 'proton'
        PrimaryIdTbl(17)%code = knuc
        PrimaryIdTbl(17)%subcode = regptcl
        PrimaryIdTbl(17)%charge = 1

        PrimaryIdTbl(18)%symb = 'p'
        PrimaryIdTbl(18)%code = knuc
        PrimaryIdTbl(18)%subcode = regptcl
        PrimaryIdTbl(18)%charge = 1

        PrimaryIdTbl(19)%symb = 'n'
        PrimaryIdTbl(19)%code = knuc
        PrimaryIdTbl(19)%subcode = kneutron
        PrimaryIdTbl(19)%charge = 0

        PrimaryIdTbl(20)%symb = 'neutron'
        PrimaryIdTbl(20)%code = knuc
        PrimaryIdTbl(20)%subcode = kneutron
        PrimaryIdTbl(20)%charge = 0

        PrimaryIdTbl(21)%symb = 'rho'
        PrimaryIdTbl(21)%code = krho
        PrimaryIdTbl(21)%subcode = 0
        PrimaryIdTbl(21)%charge = 0

        PrimaryIdTbl(22)%symb = 'omega'
        PrimaryIdTbl(22)%code = komega
        PrimaryIdTbl(22)%subcode = 0
        PrimaryIdTbl(22)%charge = 0

        PrimaryIdTbl(23)%symb = 'phi'
        PrimaryIdTbl(23)%code = kphi
        PrimaryIdTbl(23)%subcode = 0
        PrimaryIdTbl(23)%charge = 0

        PrimaryIdTbl(24)%symb = 'd+'
        PrimaryIdTbl(24)%code = kdmes
        PrimaryIdTbl(24)%subcode = regptcl
        PrimaryIdTbl(24)%charge = 1

        PrimaryIdTbl(25)%symb = 'd-'
        PrimaryIdTbl(25)%code = kdmes
        PrimaryIdTbl(25)%subcode = antip
        PrimaryIdTbl(25)%charge = -1

        PrimaryIdTbl(26)%symb = 'd0'
        PrimaryIdTbl(26)%code = kdmes
        PrimaryIdTbl(26)%subcode = regptcl
        PrimaryIdTbl(26)%charge = 0

        PrimaryIdTbl(27)%symb = 'dd~'
        PrimaryIdTbl(27)%code = kddb
        PrimaryIdTbl(27)%subcode = 0
        PrimaryIdTbl(27)%charge = 0

        PrimaryIdTbl(28)%symb = 'nn~'
        PrimaryIdTbl(28)%code = knnb
        PrimaryIdTbl(28)%subcode = 0
        PrimaryIdTbl(28)%charge = 0

        PrimaryIdTbl(29)%symb = 'neu_e'
        PrimaryIdTbl(29)%code = kneue
        PrimaryIdTbl(29)%subcode = regptcl
        PrimaryIdTbl(29)%charge = 0

        PrimaryIdTbl(30)%symb = 'neu_mu'
        PrimaryIdTbl(30)%code = kneumu
        PrimaryIdTbl(30)%subcode = regptcl
        PrimaryIdTbl(30)%charge = 0

        PrimaryIdTbl(31)%symb = 'deuteron'
        PrimaryIdTbl(31)%code = kgnuc
        PrimaryIdTbl(31)%subcode =2
        PrimaryIdTbl(31)%charge = 1

        PrimaryIdTbl(32)%symb = 'alfa'
        PrimaryIdTbl(32)%code = kalfa
        PrimaryIdTbl(32)%subcode = regptcl
        PrimaryIdTbl(32)%charge = 1

        PrimaryIdTbl(33)%symb = 'he'
        PrimaryIdTbl(33)%code = kalfa
        PrimaryIdTbl(33)%subcode = regptcl
        PrimaryIdTbl(33)%charge = 1

        PrimaryIdTbl(34)%symb = 'alpha'
        PrimaryIdTbl(34)%code = kalfa
        PrimaryIdTbl(34)%subcode = regptcl
        PrimaryIdTbl(34)%charge = 1          ! charge is internally set

        PrimaryIdTbl(35)%symb = 'l'
        PrimaryIdTbl(35)%code = klibe
        PrimaryIdTbl(35)%subcode = regptcl
        PrimaryIdTbl(35)%charge = 1

        PrimaryIdTbl(36)%symb = 'libeb'
        PrimaryIdTbl(36)%code = klibe
        PrimaryIdTbl(36)%subcode = regptcl
        PrimaryIdTbl(36)%charge = 1

        PrimaryIdTbl(37)%symb = 'm'
        PrimaryIdTbl(37)%code = kcno
        PrimaryIdTbl(37)%subcode = regptcl
        PrimaryIdTbl(37)%charge = 1

        PrimaryIdTbl(38)%symb = 'cno'
        PrimaryIdTbl(38)%code = kcno
        PrimaryIdTbl(38)%subcode = regptcl
        PrimaryIdTbl(38)%charge = 1

        PrimaryIdTbl(39)%symb = 'h'
        PrimaryIdTbl(39)%code = khvy
        PrimaryIdTbl(39)%subcode = regptcl
        PrimaryIdTbl(39)%charge = 1

        PrimaryIdTbl(40)%symb = 'namgsi'
        PrimaryIdTbl(40)%code = khvy
        PrimaryIdTbl(40)%subcode = regptcl
        PrimaryIdTbl(40)%charge = 1

        PrimaryIdTbl(41)%symb = 'vh'
        PrimaryIdTbl(41)%code = kvhvy
        PrimaryIdTbl(41)%subcode = regptcl
        PrimaryIdTbl(41)%charge = 1

        PrimaryIdTbl(42)%symb = 'sclar'
        PrimaryIdTbl(42)%code = kvhvy
        PrimaryIdTbl(42)%subcode = regptcl
        PrimaryIdTbl(42)%charge = 1

        PrimaryIdTbl(43)%symb = 'fe'
        PrimaryIdTbl(43)%code = kiron
        PrimaryIdTbl(43)%subcode = regptcl
        PrimaryIdTbl(43)%charge = 1

        PrimaryIdTbl(44)%symb = 'iron'
        PrimaryIdTbl(44)%code = kiron
        PrimaryIdTbl(44)%subcode = regptcl
        PrimaryIdTbl(44)%charge = 1

        PrimaryIdTbl(45)%symb = 'iso'
        PrimaryIdTbl(45)%code = kgnuc
        PrimaryIdTbl(45)%subcode = 3    ! to be determined later
        PrimaryIdTbl(45)%charge = 2   ! tbd

        PrimaryIdTbl(46)%symb = 'sigma0'
        PrimaryIdTbl(46)%code = ksigma
        PrimaryIdTbl(46)%subcode = -1
        PrimaryIdTbl(46)%charge = 0

        PrimaryIdTbl(47)%symb = 'sigma+'
        PrimaryIdTbl(47)%code = ksigma
        PrimaryIdTbl(47)%subcode = -1
        PrimaryIdTbl(47)%charge = 1

        PrimaryIdTbl(48)%symb = 'sigma-'
        PrimaryIdTbl(48)%code = ksigma
        PrimaryIdTbl(48)%subcode = -1
        PrimaryIdTbl(48)%charge = -1

        PrimaryIdTbl(49)%symb = 'gzai0'
        PrimaryIdTbl(49)%code = kgzai
        PrimaryIdTbl(49)%subcode = -1
        PrimaryIdTbl(49)%charge = 0

        PrimaryIdTbl(50)%symb = 'gzai+'
        PrimaryIdTbl(50)%code = kgzai
        PrimaryIdTbl(50)%subcode = -1
        PrimaryIdTbl(50)%charge = 1

        PrimaryIdTbl(51)%symb = 'bomega'
        PrimaryIdTbl(51)%code = kbomega
        PrimaryIdTbl(51)%subcode = -1
        PrimaryIdTbl(51)%charge = -1

        PrimaryIdTbl(52)%symb = 'lambda'
        PrimaryIdTbl(52)%code = klambda
        PrimaryIdTbl(52)%subcode = -1
        PrimaryIdTbl(52)%charge = 0

        PrimaryIdTbl(53)%symb = 'lambdac'
        PrimaryIdTbl(53)%code = klambdac
        PrimaryIdTbl(53)%subcode = -1
        PrimaryIdTbl(53)%charge = 1

        PrimaryIdTbl(54)%symb = 'pbar'
        PrimaryIdTbl(54)%code = knuc
        PrimaryIdTbl(54)%subcode = antip
        PrimaryIdTbl(54)%charge = -1

        
        PrimaryIdTbl(55)%symb = 'nbar'
        PrimaryIdTbl(55)%code = knuc
        PrimaryIdTbl(55)%subcode = antip
        PrimaryIdTbl(55)%charge = 0

        PrimaryIdTbl(56)%symb = 'light'
        PrimaryIdTbl(56)%code = klight
        PrimaryIdTbl(56)%subcode = 1
        PrimaryIdTbl(56)%charge = 0

        PrimaryIdTbl(57)%symb = 'Edepo'
        PrimaryIdTbl(57)%code =  kEdepo
        PrimaryIdTbl(57)%subcode = 0
        PrimaryIdTbl(57)%charge = 0

        PrimaryIdTbl(58)%symb = 'chgPath'
        PrimaryIdTbl(58)%code =  kchgPath
        PrimaryIdTbl(58)%subcode = 0
        PrimaryIdTbl(58)%charge = 1

        PrimaryIdTbl(59)%symb = 'tau+'
        PrimaryIdTbl(59)%code =  ktau
        PrimaryIdTbl(59)%subcode = antip
        PrimaryIdTbl(59)%charge = 1

        PrimaryIdTbl(60)%symb = 'tau-'
        PrimaryIdTbl(60)%code =  ktau
        PrimaryIdTbl(60)%subcode =  regptcl
        PrimaryIdTbl(60)%charge = -1

        PrimaryIdTbl(61)%symb = 'neu_tau'
        PrimaryIdTbl(61)%code =  kneutau
        PrimaryIdTbl(61)%subcode =  regptcl
        PrimaryIdTbl(61)%charge = 0

        PrimaryIdTbl(62)%symb = 'neu_e~'
        PrimaryIdTbl(62)%code = kneue
        PrimaryIdTbl(62)%subcode = antip
        PrimaryIdTbl(62)%charge = 0

        PrimaryIdTbl(63)%symb = 'neu_mu~'
        PrimaryIdTbl(63)%code = kneumu
        PrimaryIdTbl(63)%subcode = antip
        PrimaryIdTbl(63)%charge = 0

        PrimaryIdTbl(64)%symb = 'neu_tau~'
        PrimaryIdTbl(64)%code =  kneutau
        PrimaryIdTbl(64)%subcode =  antip
        PrimaryIdTbl(64)%charge = 0

        PrimaryIdTbl(65)%symb = 'pdg'
        PrimaryIdTbl(65)%code = 0       ! later fixed
        PrimaryIdTbl(65)%subcode = 0    ! //
        PrimaryIdTbl(65)%charge = 0     !

        PrimaryIdTbl(66)%symb = 'p~'
        PrimaryIdTbl(66)%code = knuc
        PrimaryIdTbl(66)%subcode = antip
        PrimaryIdTbl(66)%charge = -1

        PrimaryIdTbl(67)%symb = 'n~'
        PrimaryIdTbl(67)%code = knuc
        PrimaryIdTbl(67)%subcode = antip
        PrimaryIdTbl(67)%charge = 0
        
        
       end
!***************************************************
       subroutine csetEUnitTbl
         implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"

         ErgUnitTbl(1)%symb = 'ev'      
         ErgUnitTbl(1)%togev = 1.d-9
        
         ErgUnitTbl(2)%symb = 'kev'
         ErgUnitTbl(2)%togev = 1.d-6
        
         ErgUnitTbl(3)%symb = 'mev'
         ErgUnitTbl(3)%togev = 1.d-3
        
         ErgUnitTbl(4)%symb = 'gev'
         ErgUnitTbl(4)%togev = 1.d0
        
         ErgUnitTbl(5)%symb = 'tev'
         ErgUnitTbl(5)%togev = 1.d3
        
         ErgUnitTbl(6)%symb = 'pev'
         ErgUnitTbl(6)%togev = 1.d6
        
         ErgUnitTbl(7)%symb = 'eev'
         ErgUnitTbl(7)%togev = 1.d9
       end
