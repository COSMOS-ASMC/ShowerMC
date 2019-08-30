!c            test ksx2
!         open(13,file='c2s5001.#gd.data')
!         do 100 i=1, 50000
!            call ksx2(19., x2)
!            write(13) x2
! 100     continue
!      end
!      ***************************************************************
!      *
!      * ksx2: samples a chi-square random variable with given degree
!      *       of freedom.
!      *
!      ***************************************************************
!
!
!     /usage/ call ksx2(df, x2)
!
!  df: input.  real*4.  degree of freedom.
!  x2: output. real*4.  sampled chi-square variable
!
          subroutine ksx2(df, x2)
             call ksgmrs(df/2-1., 1., y)
             x2=2*y
          end
