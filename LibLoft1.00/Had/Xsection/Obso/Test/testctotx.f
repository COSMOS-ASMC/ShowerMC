!       sigma tot, inel, ela test
#include "ZcosmosBD.h"
      program testxs
      implicit none      
! #include "ZcosmosExt.h"
#include "Zcode.h"
#include "Zptcl.h"
!
!
      real*8  p, A, Z, Ek, xsT, XsI, XsE, xsxpt, xsxpin
      integer code, subcode, charge, i
      type(ptcl):: pj
      write(0,*) 'enter proj code, subcode, charge'
      read(*,*)  code, subcode, charge
      write(0,*)' Enter A of target '
      read(*,*) A
      Z = max(A/2.d0, 1.d0)   !dummy


      Ek = 0.01
      call cmkptc(code, subcode, charge, pj)
      write(*,*) '# Ek(GeV)  St Sin Sel  P(GeV) Sxpt Sxpin'
      write(*,'(a, 3i4)')
     *    '# pj ', pj.code, pj.subcode, pj.charge
      write(*,'(a, f8.3)') '# tg A', A

      do while ( Ek .lt. 1.d9 )
         p = sqrt( (Ek+ pj.mass)**2 - pj.mass**2)
         pj.fm.p(1) = 0.
         pj.fm.p(2) = 0.
         pj.fm.p(3) = p
         pj.fm.p(4) = Ek + pj.mass
         call ctotx(pj, A, xsT)
         call cinelx(pj, A, Z, xsI)
         call ctotx(pj, 1.0d0, xsxpt)
         call cinelx(pj, 1.0d0, 1.d0,  xsxpin)
         xsE = xsT- xsI
         write(*,'(1p,7g13.4)')  Ek, xsT, xsI, xsE, p, 
     *            xsxpt, xsxpin
         Ek = Ek*10.**0.02
      enddo

      end
