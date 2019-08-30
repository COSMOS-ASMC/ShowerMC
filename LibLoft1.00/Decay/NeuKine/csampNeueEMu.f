!       ****************************************************************
!       *
!       * csampNeueEMu:  sample energy of electron neutrino from mu decay
!       *
!       **************************tested 88.07.26***************k.k*****
!
! /usage/   call csampNeueEMu(f)
!   f: output. real*8.  sampled fractional energy.  f is the
!                       fraction  given by f=2*e'/m  where e' is the
!                       energy in muon rest system and m the muon mass.
!         
!       df f^2 (1-f)  which can be sampled by taking  the second 
!       max of 4 uniform random number.
!
          subroutine csampNeueEMu(f)
          implicit none
          real*8 f

          real*8  u1, u2, u3, u4, ux, uy
          call rndc(u1)
          call rndc(u2)
          call rndc(u3)
          call rndc(u4)
          ux = max(u1, u2)
          uy = max(u3, u4)
          if(ux .gt. uy) then
             f=max( min(u1, u2), uy)
          else
             f=max( min(u3, u4), ux)
          endif
          end
