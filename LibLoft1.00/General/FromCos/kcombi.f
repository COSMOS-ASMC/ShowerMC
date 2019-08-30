!c        test kcombi
!         program test_kcombi
!         implicit none
!         real*8 m, n, c
!         m= 10.
!         n= 10.
!         call kcombi(n, m, c)
!         write(*, *) ' n=',n, ' m=', m, ' c=',c
!         end
!          **************************************************
!          *
!          * kcombi: compute  n!/(n-m)!/m!
!          *
!          **************************************************
!
!   /usage/  call kcombi(n, m, c)
!  n:  real*8. input.
!  m:  real*8. input. must be  n>=m.
!  c:  real*8. output.  n!/(n-m)!/m!
!
!   note:  n and m may be fractional 
!
         subroutine kcombi(n, m, c)
         implicit none
         real * 8 n, m, c

         real*8 big, sq2pi
         parameter (big = 20., sq2pi =2.5066282746310)
         real*8 ep, x, z1, z2, z, tmp, kgamma

!
         ep(x)=(1.d0/288./x + 1.d0/12.)/x + 1.
!
            z1=n+1.
            z2=m+1.
            z=n-m+1.
            if(n .gt. big .and. (n-m) .gt. big .and. m .gt. big) then
                tmp=exp( (n+.5)*log( z1/z) + m*log(z/z2) + 1.
     *             -.5 *log(z2) )
                c=tmp/sq2pi * ep(z1)/ep(z2)/ep(z)
            elseif(n .gt. big .and. (n-m) .gt. big) then
                tmp=exp(m * (log(z)-1.) + (n+.5)*log(z1/z) )
                c=tmp/kgamma(z2)*ep(z1)/ep(z)
            elseif(n .gt. big .and. m .gt. big) then
                tmp=exp( (n+.5)*log(z1) - (m+.5)*log(z2) + m-n)
                c=tmp*ep(z1)/kgamma(z)/ep(z2)
            else
                c=kgamma(z1)/kgamma(z2)/kgamma(z)
            endif
        end
