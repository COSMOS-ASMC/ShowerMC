!c     test ksgamd
!      include 'kgauss.f'
!      include 'kgamma.f'
!      include 'rnd.f'
!--------------------------------
!     implicit none
!     real*8 alfa, am, xm, x
!     integer i
!
!        alfa=15.
!        am=105.
!        xm=alfa/am/(am-alfa)
!        do  i=1, 10000
!c            activate one of the following calls.
!c          call ksgmim(3, 300.d0, x)
!c          call ksgmis(2, 10.d0, x)
!c          call ksgmrm(3.5d0, 100.d0, x)
!c          call ksgmrs(3.5d0,  20.d0, x)
!c          call ksgmrs(0.2d0, 30.d0, x)
!c          call ksgmrs(-0.45d0, 300.d0, x)
!c          call ksgmrs(-0.95d0, 20.d0, x)
!c          call ksgmrm(1.00d0, 2.d0,  x)
!           write(*, *)  sngl(x)
!        enddo  
!      end
!    *****************************************************************
!    * sampling for the gamma distribution: (note power is not s-1)
!    *
!    *  (x/a)**s exp(-x/a)/gamma(s+1)d(x/a).
!    *
!    *     (s>-1.0; mean is a*(s+1)).
!    *
!    * ksgmim: sampling for an integer s with a given mean
!    * ksgmis: sampling for an integer s with a given scale parameter a
!    * ksgmrm: sampling for a real s with a given mean
!    * ksgmrs: sampling for a real s with a given scale parameter a
!    *
!    ****************************************************************
!
!    samples a random variable from gamma distribution, x**(s) exp(-x/a)
!
!   /usage/
!    call ksgmim(n, av, x)
!            n: input. integer>= 0.   s=n.
!           av: input.  real*8> 0.   the average of the distribution
!            x: output. real*8.  sampled random variable.
!
!    call ksgmis(n, a, x)
!            n: input. integer>=0.   s=n
!            a: input. real*8.    a in exp(-x/a)
!                      the mean is given by a*(n+1)
!            x: output. real*8.  sampled random variable.
!
!    call ksgmrm(s, av, x)
!            s: input. real*8 >-1.0
!           av: input. real*8.  >0. the average of the distribution.
!            x: output. real*8. a sampled radonm varible.
!                      if s is close to -1.0, x < 1.d-30 may not
!                      be very accurate.
!
!    call ksgmrs(s, a, x)
!            s: input.  real*8 > -1.
!            a: input.  real*8 > 0. a in exp(-x/a)
!                      the mean is given by a*(s+1)
!            x: output. real*8.  sampled random variable.
!                      if s is  close to -1.0, x < 1.d-30 may not
!                      be very accurate.
!
!         for s=3.5, 1653 msec is needed for calling 50000 times.
!         (by m780 35mips machine)
!
          subroutine ksgmrs(s, a, x)
!           implicit none
            real*8 s,  a, x
!
            real*8 alfa, u, r
            integer n
            logical more
!
              if(s .le. -1.) then
                    write(*,*) ' input value to ksgmrs. s=',
     *              s,' invalid'
                    stop 9999
              elseif(s .lt. 0.) then
                    call ksgamn(s, x)
              else
!                     sample a varible from x**s exp(-x)dx
!                      decompose s=n+alfa,  1.0>alfa>=0.
                    n=s
                    alfa=s-n
                    if(alfa .eq. 0.) then
!                        use integer case
                        call ksgmis(n, 1.d0,x)
                    else
!                        importance sampling with weight function
!                        rho= (1-alfa)*x**n exp(-x)/gamma(n+1) +
!                        alfa*x**(n+1) exp(-x)/gamma(n+2)
!                        =x**n exp(-x)/gamma(n+1) (1-alfa+alfa*x/(n+1))
!                        then,
!                        x**(n+alfa)exp(-x)/rho=
!                        x**alfa/(1-alfa+alfa*x/(n+1))
!                        which has a max value, (n+1)**alfa at x=n+1.
!                        determine if x**n exp(-x) or x**(n+1) exp(-x)
                        more=.true.
                        do 100 while (more)
!                            average trial # is 1.03 with max of 4 (for
!                            s=3.5, 50000 cases).
                            call rndc(u)
                            if(u .le. alfa) then
                                 call ksgmis(n+1,1.d0, x)
                            else
                                 call ksgmis(n, 1.d0, x)
                            endif
!                             rejection func.value
                            r= (x/(n+1))**alfa/(1.-alfa + alfa*x/(n+1))
                            call rndc(u)
                            more=u .gt. r
  100                   continue
                   endif
              endif
              x=a*x
          end
!
          subroutine ksgmrm(s, av, x)
!         implicit none
          real*8 s, av, x
!
          real*8 a
!
             if(s .le. -1.) then
                    write(*,*) ' input value to ksgmrs. s=',
     *              s,' invalid'
                    stop 9999
             endif
             a=av/(s+1.)
             call ksgmrs(s, a, x)
          end
!
          subroutine ksgmis(n, a, x)
!         implicit none
          real*8 a, x
          integer n
!
          real*8 u, ui
          integer i
          logical more
!
             if(n .le. -1) then
                  write(*,*) ' input value to ksgmis. n=',
     *            n,' invalid'
                  stop 9999
             endif
!
             more=.true.
             do 100 while (more)
                 call rndc(u)
                 do 50  i=1, n
                     call rndc(ui)
                     u=u*ui
   50            continue
                 more=u .le. 0.
  100        continue
             x=-a*log(u)
           end
!
           subroutine ksgmim(n, av, x)
!          implicit none
           real*8 av, x
           integer n
!
           real*8 a
!
             if(n .le. -1) then
                  write(*,*) ' input value to ksgmim. n=',n,' invalid'
                  stop 9999
             endif
!
             a=av/(n+1)
             call ksgmis(n, a, x)
           end
!            for   -1.0<s<0.
           subroutine ksgamn(s, x)
!          implicit none
           real*8 s, x
!
!
           real*8 pi/3.14159265/, rt2i/0.70710678/, eps/0.1d0/,
     *           err/1.d-3/,   small/1.d-30/
           real*8 u, r, us, br, kgamma, tmp, xold, acc
           logical more, more2
!
             if(s .ge. -0.5) then
!    importance sampling is completely possible by using weight
!    function of
!    rho=-2s* x**(-1/2)exp(-x)/gamma(1/2)+(2s+1)exp(-x).
!    then, x**s exp(-x)/rho=
!    x**s/( -2s x**(-1/2)/gamma(1/2) + (2s+1)).
!    this takes the
!    maximum, (1/gamma(1/2))**2s=pi**(-s) at x=(1/gamma(1/2))**2=1/pi
!    sampling for exp(-x)/root(x)dx  is done by taking the square of
!    a gaussian random varible with mean 0 and variance 1/root(2).
!
                 more=.true.
                 do 100 while ( more )
                    call rndc(u)
                    if(u .le. (2*s+1.)) then
                        call rndc(u)
                        x=-log(u)
                    else
                        call kgauss(0.d0, rt2i, x)
                        x=x*x
                    endif
!                         rejection function value
                    r=(pi*x)**(s+0.5)/(-2*s + (2*s+1.)*sqrt(pi*x))
                    call rndc(u)
                    more=u .gt. r
  100           continue
             else
!                 complete rejection method not known for me(ptcl data
!                 show some, but seems wrong).
!                 introduce a cut at a small x(=eps).
!                 divide the samling region into 2: x<eps and x> eps
!                 at x>eps
!                  use weight function exp(-x)/root(x)
!                  x**s exp(-x)/rho= x**(s+1/2)
!                at x<eps
!                  solve the usual equation approxmately. for small x,
!                  int(0;x)(x**s exp(-x)/gamma(s+1))=
!                  x**(s+1)/gamma(s+2)(1-(s+1)/(s+2)x+(s+1)/(s+3)/2x**2
!                 then,  the relative area of x< eps
                br=eps**(s+1.)/kgamma(s+2.) * (1.- (s+1.)/(s+2.)*eps
     *             + (s+1.)/(s+3.)/2*eps*eps  )
                call rndc(u)
                if(u .lt. br) then
!                    take value x< eps
                   call rndc(u)
                   us=u**(1./(s+1.))*eps
!                     initial guess at x<eps
                   x=us
                   tmp=
     *             (1.- (s+1.)/(s+2.)*eps+(s+1.)/(s+3)/2*eps*eps)
     *             **(1./(s+1.))
                   more=.true.
                   do 130  while( more )
!                         average trial is < 2
                       xold=x
                       x=us*tmp/
     *                 (1.-(s+1.)/(s+2.)*xold+(s+1.)/(s+3)/2*xold*xold)
     *                 ** (1./(s+1.))
                       if(x .lt. small) then
                          acc=0.
                       else
                          acc=abs(x/xold-1.)
                       endif
                       more=acc .gt. err
  130               continue
                else
                   more=.true.
                   do 200  while ( more )
                       more2=.true.
                       do 150 while ( more2 )
                           call kgauss(0.d0, rt2i, x)
                           x=x*x
                           more2= x .lt. eps
  150                  continue
!                             rejection function value
                       r=(x/eps)**(s+0.5)
                       call rndc(u)
                       more=u .gt. r
  200              continue
               endif
            endif
          end



