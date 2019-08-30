      subroutine cgetMFP(sigma, rhoNabyA, mfp)
      implicit none
      real*8 sigma  ! input sum r_i sigma_i. r_i is the normalzied relatvie numer
                    ! of target elments.  sigma_i the cross-section of i-th element
                    !  (  cm^2)
      real*8 rhoNabyA  !input. rho*Na/<A>. in /cm^3.
                       !  <A>= sum of r_i A_i 
      real*8 mfp       ! output.  mean free path of collision. in cm.

      if(sigma .gt. 0.0d0) then
         mfp = 1.d0/(rhoNabyA*sigma)
      else
         mfp = 1.d35
      endif
      end

      subroutine cgetsigmai(A, Elab, code, sigmai)
      end

 
