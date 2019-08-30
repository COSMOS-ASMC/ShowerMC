c  ifort test.f  -L$LIBLOFT/lib/MacIFC -lloft
c This is to test the sampling scheme for prob. density function 
c  consisting of  F(x)dx= ( f(x) + g(x) )
c     where F(x) and f(x) are all >=0 but g(x) is <= 0
c     
c  1)     sample a random variable x from f(x)
c  2)     get eff from  eff* f(x)  = F(x) = f(x) + g(x)
c      i.e, eff = 1+ g(x)/f(x) < 1
c  3)     sample a random number u uniform  in  (0,1) and
c         accept x if u<eff.
c
c
      implicit none
      real(8):: u, x
      real(8):: norm, A, B, eff
      integer i
c
c        f(x) dx = [exp(-x) +x(x-2)*0.25 ]dx
c                    f(x)   + g(x)
c
c       integral exp(-x) dx 0 to 2 = -exp(-x) (0 to 2) = 1- exp(-2)
c
      norm = 1.d0 - exp(-2.d0)
c 
c       sample x from f(x)dx
c        norm * u = (1-exp(-x));  1- norm *u=exp(-x)
c
      do i = 1, 1000000
         call rndc(u)
         x = - log(1.d0- norm*u)
c                 
         A = exp(-x)
         B = x*(x-2)*0.25
c      C = A + B 
c          A * eff = C =A + B ; eff = 1+ B/A
         eff = 1. + B/A
         call rndc(u)
         if( u < eff ) then
            write(*,*) x
         endif
      enddo
      end

      
      
