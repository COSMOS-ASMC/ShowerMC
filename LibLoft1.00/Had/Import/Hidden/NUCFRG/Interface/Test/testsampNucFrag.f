      program main
      use NucMassFormula
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
      type(ptcl)::pj
      integer:: ita, itz, Nsp, A, Z
      integer,parameter::npmax=500
      type(ptcl)::x(npmax)
      
      integer:: i, np, j

      A = 103
      Z = 45
      ita = 12
      itz = 6
      call cmkptc(kgnuc, A, Z, pj)
      pj.mass = cNucMass(A, Z)/1000.
      pj.fm.p(1:2) = 0.
      pj.fm.p(3) = 5*A
      pj.fm.p(4) =sqrt( pj.fm.p(3)**2 +  pj.mass**2 )
      Nsp = 60
      do i = 1, 100000
         call csampNucFrag(pj, ita, itz, Nsp, x, np)      
         do j = 1, np 
            write(*,'(i7, 2i4,3i5, 1p, 5g14.5)' ) i, j, np,
     *     x(j).code, x(j).subcode, x(j).charge,
     *     x(j).fm.p(1:4), x(j).fm.p(4)-x(j).mass
         enddo
      enddo
      end
