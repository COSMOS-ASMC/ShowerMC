#include "BlockData/cblkGene.h"
#include "Zphitsblk.h"                                                                            
      implicit none
#include "Zcode.h"
#include "Zmass.h"
#include "Zptcl.h"
#include "ZcosmosExt.h"
      type(ptcl):: pj
      real(8)::A, xs, Z
      real*8 p, Ek, mass, Ekpn
      integer code, subc, charge
      write(0,*) " Enter pj's:  A Z"
      read(*,*)  subc, charge
      write(0,*) "Enter target A, Z"
      read(*,*)  A, Z

      code = kgnuc
      call cmkptc(code, subc, charge, pj)

      Ekpn = 0.01   ! GeV
      
      write(*,*)
     *   '# pj A,Z=', subc, charge, ' target A,Z',int(A), int(Z)
      write(*,*)  '# Ekt,  xs, p, Ek/n  p/n (E in GeV)'
      do while( Ekpn .lt. 1.e7 )
         Ek = Ekpn*subc
         mass = pj.mass
         p = sqrt( (Ek+ mass)**2 - mass**2)
         pj.fm.p(1:2) = 0.
         pj.fm.p(3) = p
         pj.fm.p(4) = Ek +mass
         call cAAXsec2(pj, A, Z, xs)
         write(*,'(1p,5g12.3)')  Ek, xs, p, Ekpn, p/subc
         Ekpn = Ekpn *10.**0.1
      enddo
      end
