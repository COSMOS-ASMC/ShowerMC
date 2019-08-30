!      compute  v K2(2zeta) which is the sencod term of
!      F in Brainerd and Petrosian's Eq(5). v is fractional energy
!      of emitted gamma ray.
!
      real*8 function cmBremF2(up,  v)
      implicit none
      real*8 up  !  input. upsilon
      real*8 v   !  input. faractional energy.
      real*8 zeta2, ck23
      
      zeta2 = v/(1.d0-v) /3.d0/up * 2.d0
      cmBremF2 = v * ck23(zeta2)*zeta2
      end
