!      real*8 xai,  x, ans
!      read(*,*) xai
!      x = 0.001d0
!      do while  (x .lt. 1.d0)
!c             devide by xai: xsection is to be divided by Eg/m
!c          but if B is const, xai ~ Eg/m so, to see the
!c        cross-section within a moderate range, we use xai division
!c         
!         ans = cmPairSpec(xai, x)/xai
!         write(*, *) sngl(x), sngl(ans)
!         x = x + 0.001d0
!      enddo
!      end
      real*8 function cmPairSpec(xai, x)
      implicit none
!         This computes electron energy distribution function
!      The value dose not contain the coefficient which is
!      SyncCoef/(3Piroot(3))/Eg -->unit is  /m
      real*8 xai  ! input. Earger's (1966) xai = Eg/m * B/Bc/2
      real*8 x    ! input. Fractional electron energy. Ee/Eg
!
      real*8 y, v, ck23

      if(x .eq. 1. .or. x .eq. 0.) then
         cmPairSpec = 0.
      else
         v = 1. - 2*x
         y = 4./(3.*xai) /(1.-v**2)
         cmPairSpec = (9.- v*v)/(1.-v*v) * ck23(y)
      endif
      end

