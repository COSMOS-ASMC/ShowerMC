!               test kbinom
!      implicit none
!      real*8  p
!      integer n0, i, n
!       p=.45
!       n0=30
!       do  i=1, 10000
!           call kbinom(p, n0,  n)
!           write(*, *) n
!       enddo  
!       end
!     ************************************************************
!     *
!     * kbinom: samples binomial random #
!     *
!     ************************ tested 88.08.26 ***********k.k*****
!
!      prob(nb)=comb(n,nb)*p**nb * (1.-p)**(n-nb)
!
! /usage/ call kbinom(p, n, nb)
!      p: input.  (0, 1)
!      n: input.  >0 integer
!     nb: output. sampled #
!
      subroutine kbinom(p, n, nb)
      implicit none
      real*8 p
      integer n, nb

      real*8 av, s, x, u
      integer i
!
          nb=0
          if(n .gt. 20 .and. p .gt. .25 .and. p .lt. .75) then
              av=p*n
              s=sqrt(p*n*(1.-p))
!              *** until loop*** 
              do while (.true.)
                 call kgauss(av, s, x)
                 nb=x+.5
                 if(nb .ge. 0 .and. nb .le. n) goto 100
              enddo
  100         continue
          else
               do   i=1, n
                  call rndc(u)
                  if(u .le. p) then
                      nb=nb+1
                  endif
               enddo
          endif
       end
