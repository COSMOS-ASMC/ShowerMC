!       ****************************************************************
!       *
!       * csNeumuEMu:  sample energy of muon neutrino from mu decay
!       *
!       **************************tested 88.07.26***************k.k*****
!
! /usage/      call csNeumuEMu(f)
!   f: output. real*8.  sampled fractional energy.  f is the
!                       fraction  given by f=2*e'/m  where e' is the
!                       energy in muon rest system and m the muon mass.
!          The angle integrated energy distribution is 
!         f^2 (3-2f)df =(2f^2(1-f) + f^2)df
!                         1:           2

        subroutine csNeumuEMu(f)
          implicit none
          real*8 f
          real*8 u, u1, u2, u3
          call rndc(u)
          if(u .lt. 0.333333) then
!                 first term
             call csampNeueEMu(f)
          else
!                 second term. take max of 3 u's.
             call rndc(u1)
             call rndc(u2)
             call rndc(u3)
             f= max(u1, u2, u3)
          endif
       end

