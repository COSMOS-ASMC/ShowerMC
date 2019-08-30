!c            test ksbeta
!         open(13,file='c2s5001.#gd.data')
!         a=1.5
!         b=2.5
!         call timei('             ', '               ', 40)
!         do 100 i=1, 50000
!            call ksbeta(a, b, x)
!            write(13) x
! 100     continue
!         call timeg(i, j, -1)
!      end
!      ***************************************************************
!      *
!      * ksbeta: samples a beta-function random variable
!      *
!      ***************************************************************
!
!
!     /usage/ call ksbeta(a, b, x)
!     x is sampled from the beta function,
!        f(x)dx=x**(a-1) * (1.-x)**(b-1) dx
!   a: input. real*4 >0.
!   b: input. real*4 >0.
! x: output. real*4.  sampled variable.
!
!   method:  let n1=2a, n2=2b, and,  x1 and x2 be chi-square random
!            variables with degree of freedom n1 and n2, resp.
!            then, x=x1/(x1+x2) follows beta distribution.
!
!  to sample 50000 variables with a=1.5, b=2.5, 3.4sec (with
!  facom m780, 35 mips) is needed.
!
!     kbetar is sometimes faster than this for a class a,b.
!
          subroutine ksbeta(a, b, x)
              df1=2*a
              df2=2*b
              call ksx2(df1,x1)
              call ksx2(df2,x2)
              x= x1/(x1+x2)
          end
