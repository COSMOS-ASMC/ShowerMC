#include "ZcosmosBD.h"
      implicit none
#include "Zptcl.h"
!#include "Zcode.h"
!#include "Zmass.h"
      real*8  xs, p, xspp, xspA, A
      integer iat
      type(ptcl):: pj

!     
      real(8):: Ek
      integer code, subcode, charge, i

      write(0,*)
     * ' Enter proj. code, subcode, charge  target A'
      read(*,*)  code, subcode, charge,iat

      Ek = 50.0
      A = iat
      call ciniQGS
      call cmkptc(code, subcode, charge, pj)
      write(*,'(a)')
!     *    '# Ek(GeV) Xs(targetA) Xs(pp) Xs(targetAbyCos) mom'
     *    '# Ek(GeV) Xs(targetA) mom'

      do while (Ek .lt. 1.e11)
         p =sqrt(  (Ek+pj.mass)**2 - pj.mass**2 )
         pj.fm.p(1) = 0.
         pj.fm.p(2) = 0.
         pj.fm.p(3) = p
         pj.fm.p(4) = Ek + pj.mass
!         Ek =  pj.fm.p(4) - pj.mass
         call cxsecQGS( pj, iat,   xs )
!         call cxsecQGS( pj, 1,  xspp)
!         call cxp2xAXsec(A, xspp,  xspA)
!         write(*,'(1p,5g12.4)')  Ek, xs, xspp, xspA, p
         write(*,'(1p,3g12.4)')  Ek, xs, p
         Ek = Ek*10.0**0.1
      enddo
      end
