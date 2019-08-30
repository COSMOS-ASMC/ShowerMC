!      include 'kgauss.f'
!      include 'rnd.f'
!c           test kpoisn
!      do i=1, 100000
!         call kpoisn(10.0d0, n)
!         write(*,*)  n
!      enddo
!      end
!c    ****************************************************************
!    *                                                              *
!    * kpoisn: Samples a poisson random variable                    *
!    *                                                              *
!    ********************** tested 86.12.18 *********************k.k*
!
!  /usage/  call kpoisn(am,  n)
!     am:  real*8. Input. The mean of the poisson distribution.
!      n:  integer output. integer.  A sampled variable.
!
!   If am >= 20, gaussian approximation is used.
!
!           subroutine needed.  rndc, kgauss.


      subroutine kpoisn(am, n)
      implicit none
      real*8 am
      integer n
!
      real*8 amsv/-1.98765d37/, ammx/20./
      real*8 avsv, psv, q, u, sqam, x
      save amsv, avsv, psv, sqam
!
      logical more
!
      if(am .lt. ammx) then
          if(amsv .ne. am) then
              avsv=am
              psv=exp(-am)
          endif
          n=-1
          q=1.
          more=.true.
          do while ( more )
              n=n+1
              call rndc(u)
              q=q*u
              more= q .ge. psv
          enddo
       else
          if(amsv .ne. am) then
              amsv=am
              sqam=sqrt(am)
          endif
          more=.true.
          do while ( more )
              call kgauss(am, sqam, x)
              n=x+.5
              more=n .lt. 0
          enddo    
       endif
       end
